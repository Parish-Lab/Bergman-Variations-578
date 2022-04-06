"""Generate & place probes for NICS computation."""

__author__ = "Dominic A. Sirianni, PhD"
__email__ = "sirianni.dom@gmail.com"
__license__ = "BSD-3-Clause"
__copyright__ = "(c) 2020-2022 D. A. Sirianni, PhD"
__date__ = "2022-04-06"

import qcelemental as qcel
import numpy as np
import argparse
import tempfile
import networkx as nx
import warnings
import inspect
import glob
import os, sys

from copy import deepcopy

# ==> Typing <==
from typing import List, Dict, Tuple, Optional, Iterable, Callable, Any, Union, NewType
import numpy.typing as npt

Molecule = qcel.models.Molecule

# ==> Module Functions <==

def xyz2mol(filename: str, dtype: str = 'xyz', **kwargs) -> Molecule:
    """Reads molecular geometry & info from XYZ file and creates QCElemental molecule.

    Lightweight wrapper around `qcel.models.Molecule.from_file()` for user convenience.

    Parameters
    ----------
    filename : str
        Name of the file to be read. 
    dtype : {'xyz', 'xyz+'}
        Format for XYZ file. `'xyz+'` denotes "enhanced" XYZ format; see Notes below.
    **kwargs : Dict, optional
        Any other parameters to pass to Molecule constructor. See qcel.models.Molecule
        documentation for full list of accepted kwargs.

    Returns
    -------
    mol : Molecule
        Molecule object for molecule specified in file `filename`.

    Notes
    -----
    1. Two file formats are accepted:

        xyz - Strict XYZ format
        -----------------------

            String Layout
            -------------
            <number of atoms>
            comment line
            <element_symbol or atomic_number> <x> <y> <z>
            ...
            <element_symbol or atomic_number> <x> <y> <z>

            QM Domain
            ---------
            Specifiable: geom, elem/elez (element identity)
            Inaccessible: mass, real (vs. ghost), elbl (user label), name, units (assumed [A]),
                          input_units_to_au, fix_com/orientation/symmetry, fragmentation,
                          molecular_charge, molecular_multiplicity

            Notes
            -----
            * <number of atoms> is pattern-matched but ignored.
            * Molecular charge assumed to be 0 (neutral)
            * Molecular multiplicity assumed to be 1 (singlet)

        xyz+ - Enhanced XYZ format
        --------------------------

            String Layout
            -------------
            <number of atoms> [<bohr|au|ang>]
            [<molecular_charge> <molecular_multiplicity>] comment line
            <psi4_nucleus_spec> <x> <y> <z>
            ...
            <psi4_nucleus_spec> <x> <y> <z>

            QM Domain
            ---------
            Specifiable: geom, elem/elez (element identity), mass, real (vs. ghost), elbl (user label),
                         units (defaults [A]), molecular_charge, molecular_multiplicity
            Inaccessible: name, input_units_to_au, fix_com/orientation/symmetry, fragmentation

            Notes
            -----
            <number of atoms> is pattern-matched but ignored.

        These format descriptions are taken from MolSSI's QCElemental project
        (https://github.com/MolSSI/QCElemental) at
        QCElemental/qcelemental/molparse/from_string.py (l. 95--136) and are reproduced here
        only for sake of convenience. All credit given to the original authors.

    2. If xyz format is specified through `dtype` argument, then the molecule is
        assumed to be a neutral singlet (charge=0, multiplicity=1). However, these
        parameters may be toggled by explicitly passing `molecular_charge` and
        `molecular_multiplicity` kwargs.
    """
    # Set default chg/mult specs if plain XYZ detected
    options = {}
    if dtype == 'xyz':
        options['molecular_charge'] = 0
        options['molecular_multiplicity'] = 1

    # Set molecule name
    filebase = filename.split('/')[-1].split('.')[0]
    options['name'] = filebase
    
    # Override defaults if chg/mult explicitly passed above
    options.update(kwargs)

    # Create molecule from XYZ string repr
    mol = qcel.models.Molecule.from_file(filename, dtype=dtype, **options)

    return mol

def find_rings(mol: Molecule) -> Union[npt.ArrayLike, None]:
    """Finds unique, simple, elementary rings in a molecule.

    Parameters
    ----------
    mol : Molecule
        Molecule within which to find rings

    Returns
    -------
    rings : Union[npt.ArrayLike, None]
        None if no rings found, List[List[int]] if found.
    """
    # Get connectivity map from BFS search based on vdW radii
    connectivity = qcel.molutil.connectivity.guess_connectivity(mol.symbols, mol.geometry)

    # Build molecular graph from connectivity (edge list in NX terminology)
    G = nx.Graph(connectivity)
    
    # Cast half-graph G to fully connected digraph
    D = G.to_directed()

    # Find simple cycles of G
    cycles = nx.simple_cycles(D)

    # Prune out self & bond cycles. Cast remaining cycles to tuple so it plays nice with set
    pruned = [tuple(c) for c in cycles if len(c) >= 3]

    # Do we have nontrivial cycles?
    if len(pruned) == 0:
        return None
    
    # Find unique cycles
    unique = set()
    if len(pruned) == 2:
        if cycle_sameness(pruned[0], pruned[1]):
            unique.add(pruned[0])

    elif len(pruned) > 2:
        for i in range(len(pruned)):
            for j in range(i+1,len(pruned)):
                if cycle_sameness(pruned[i], pruned[j]):
                    unique.add(pruned[i])

    # If we only have one unique cycle, return it
    if len(unique) == 1:
        return list(unique)

    # If we have more than one unique cycle, only keep elementary ones (no subcycles)
    # Test cycles for presence of _chords_, i.e., edges b/t vertices _not_ contained in cycle
    chorded = set()
    for i in unique:
        i_set = set(i)
        for j in unique:
            j_set = set(j)
            # Test using proper subset comparison
            if i_set < j_set:
                chorded.add(j)
            elif j_set < i_set:
                chorded.add(i)

    # Take set difference of unique - chorded for elementary cycles
    elementary = unique - chorded

    # Typecasting elementary from Set[Tuple] -> List[List]
    rings = list(map(list, elementary))

    return rings

def cycle_sameness(first: list, second: list) -> bool:
    """Evaluates the sameness of two graph cycles.

    Parameters
    ----------
    first : list    
        First graph cycle against which to compare.
    second : list
        Second graph cycle to compare

    Returns
    -------
    same : bool
        True if the cycles are the same up to circular permutation, False if not
    """
    # If the cycles aren't the same length, they're not the same
    if len(first) != len(second):
        return False

    # If they don't contain the same elements, they're not the same
    if len(set(first) - set(second)) != 0: 
        return False

    # Now test for circular similarity using string subsetting
    str1 = ' '.join(map(str, first))
    dstr2 = ' '.join(map(str, second*2))
    rdstr2 = ' '.join(map(str, (second*2)[::-1]))

    if str1 in dstr2:
        return True
    # Order might be reversed, check for that
    elif str1 in rdstr2:
        return True
    else:
        return False

def principal_inertial_axis(mol: Molecule,
                            subset: Union[npt.ArrayLike, None] = None,
                            masswt: bool = True) -> Tuple[npt.NDArray, npt.NDArray]:
    """Gets inertial axis & centroid for molecule `mol`.

    Parameters
    ----------
    mol : Molecule
        Molecule for which to evaluate the inertial axis
    subset : npt.ArrayLike or NoneType, default=None
        List of atom indices for which to get inertial axis.
    masswt : bool, default=True
        Mass weight coordinates before getting inertial axis?

    Returns
    -------
    centroid : npt.NDArray
        Centroid of molecular subset; depends on masswt. See Note 1 below.
    axis : npt.NDArray
        Unit vector in R^3 defining the inertial axis with origin at molecule's centroid

    Notes
    -----
    1. Centroid location varies depending on the choice of `masswt`. If `masswt=True`,
        then the centroid will coincide with the molecular center of mass (COM). Otherwise,
        the centroid will simply be the geometric mean of specified atomic positions.
    """

    # Build weighting vector
    if masswt:
        weights = np.asarray(mol.masses)
    else:
        weights = np.ones(mol.masses.shape)
    
    # Subset of the molecule?
    if subset is not None:
        coords = np.take(mol.geometry, subset, axis=0)
        weights = np.take(weights, subset)

    else:
        coords = mol.geometry.copy()

    # Translate molecule (or subset) to centroid?
    centroid = np.average(coords, axis=0, weights=weights)

    # Build inertial tensor
    tensor = qcel.models.Molecule._inertial_tensor(coords, weight=weights)

    # Solve inertial system & return principal axis
    _, ev = np.linalg.eigh(tensor)

    return centroid, ev[:,-1]

def place_probe_along_axis(centroid: npt.ArrayLike, axis: npt.ArrayLike, dist: float
                          ) -> npt.NDArray[float]:
    """Receives array_like defining centroid & axis and places probe at location `loc`.

    Parameters
    ----------
    centroid : npt.ArrayLike
        Point to use as origin in NICS probe placement
    axis : npt.ArrayLike
        Axis along which to place probes
    dist : float
        Distance from `centroid` along `axis` to place probe

    Returns
    -------
    probe : npt.ArrayLike
        (1,3) point in Cartesian space (<x, y, z>) where probe is to be placed.
    """
    # From `centroid`, go `dist` along `axis` to generate probe location
    probe = centroid + (axis * dist)

    return probe

def augment(mol: Molecule, centroid: npt.ArrayLike, axis: npt.ArrayLike,
            recipe: str = None, loc: npt.ArrayLike = None) -> Molecule:
    """Augments XYZ geom with NICS probes placed along `axis` for Molecule `mol` centered
        at `centroid`, according either to `recipe` or at provided locations `loc`.

    Parameters
    ----------
    mol : Molecule
        Complete Molecule (not a subset) to augment with NICS probes
    centroid : npt.ArrayLike
        Ring centroid to use as origin for inertial system & NICS probe placement
    axis : npt.ArrayLike
        Axis along which to place NICS probes.
    recipe : {None, 'IsoSimple', 'IsoZScan'}
        Recipe by which to place NICS probes
        Options:
            - `'IsoSimple'` : NICS(0) + NICS(-1) + NICS(+1)
            - `'IsoZScan'` : Scan along axis for isotropic shielding at [1.0, 5.0, 0.1]
    loc : npt.ArrayLike or NoneType, default=None
        Explicit locations at which to place NICS probes. Does not have to be along `axis`.
    write : bool, default=False
        Write augmented molecule to XYZ file?
    filename : str, default='nics.xyz'
        Name of XYZ file to write augmented molecule to
    probe_symbol : {'X', 'GH', 'He'}
        Symbol to define NICS probe(s) in XYZ file. See note 1 below for details on options.
    dtype : {'xyz+', 'xyz'}

    Returns
    -------
    augmented : Molecule
        Molecule object augmented with NICS probes.

    Notes
    -----
    1. The appropriate probe symbol will depend on the downstream use of the XYZ file.
    Possible applications & recommended symbols:
        - Internal/Psi4/QCElemental format: 'X'
        - Q-Chem format: 'GH'
        - Visualization with IQMol: 'He'
    2. This version of the function works by generating an intermediate molstr representation,
    instead of augmenting a schema & remaking a mol directly. Not my circus, not my monkeys.
    """
    # Sanity check: are there either a recipe and/or explicit locations passed?
    if recipe is None and loc is None:
        raise Exception("Neither `recipe` nor `loc` requested! Please pass probe(s).")

    # Cast `loc` as list so we can add points from the recipe
    if loc is not None:
        loc = list(loc)
    else:
        loc = []

    # ==> Recipe <==
    # Interpret NICS recipe
    if recipe is not None:
        axis_probes = []
        if recipe.lower() == 'isosimple':
            axis_probes += [-1.0, 0.0, 1.0]
        elif recipe.lower() == 'isozscan':
            axis_probes += list(np.arange(-5.0, 5.0, 0.1))

        # Place probes according to recipe
        for dist in axis_probes:
            # Convert dist from Ang->Bohr
            dist *= qcel.constants.conversion_factor("angstrom", "bohr")
            loc.append(list(place_probe_along_axis(centroid, axis, dist)))

    # ==> Augment `mol` with probes <==
    # Cast molecule to schema so it's mutable
    molrec = deepcopy(mol.dict())

    # Add probes to geom & symbols list
    geom_list = list(molrec['geometry'])
    symbols_list = list(molrec['symbols'])

    for probe in loc:
        geom_list.append(probe)
        # Internal: Use `X` for probe
        symbols_list.append('X')

    # Make molrec dict for feeding to_string()
    molrec['geom'] = np.asarray(geom_list)
    molrec['symbols'] = np.asarray(symbols_list)
    # to_string() needs this extra stuff
    molrec['elem'] = np.asarray(symbols_list)
    molrec['elea'] = [-1] * len(symbols_list)
    molrec['elez'] = [-1] * len(symbols_list)
    molrec['mass'] = [-1] * len(symbols_list)
    molrec['elbl'] = [-1] * len(symbols_list)
    molrec['real'] = [True] * len(symbols_list)
    molrec['units'] = 'Bohr'

    # Build new molecule instance from temporary file.
    molstr = qcel.molparse.to_string(molrec, dtype='xyz', units='Bohr')
    with tempfile.NamedTemporaryFile('w+t') as temp:
        temp.write(molstr)
        temp.seek(0)
        augmented = qcel.models.molecule.Molecule.from_file(temp.name, dtype='xyz+')

    return augmented

def sanitize_augmented(augmented: Molecule, probe_symbol: str = 'He', dtype='xyz') -> str:
    """Sanitize a molecule augmented with NICS probes to format defined by `probe_symbol`.

    Parameters
    ----------
    augmented : Molecule
        Molecule augmented with NICS probes for which to write to file.
    probe_symbol : {'He', 'GH', 'Bq', 'X'}
        Symbol to define NICS probe(s) in molstr
    dtype : {'xyz', 'qchem'}
        Format for molstr.

    Returns
    -------
    molstr : str
        Sanitized string repr of molecule 

    Notes
    -----
    1. The appropriate probe symbol will depend on the downstream use of the XYZ file.
    Possible applications & recommended symbols:
        - Visualization with IQMol: 'He'
        - Q-Chem format: 'GH'
        - Gaussian16/GaussView format: 'Bq'
        - Internal/Psi4/QCElemental format: 'X'
    """
    # Cast molecule to schema so it's mutable & we can change probe symbols
    molrec = deepcopy(augmented.dict())

    # Add probes to geom & symbols list
    geom_list = list(molrec['geometry'])
    symbols_list = list(molrec['symbols'])

    # Make molrec dict for passing -> to_string()
    molrec['geom'] = np.asarray(geom_list)
    # to_string() needs this extra bullshit
    molrec['elea'] = [-1] * len(symbols_list)
    molrec['elez'] = [-1] * len(symbols_list)
    molrec['mass'] = [-1] * len(symbols_list)
    molrec['elbl'] = [-1] * len(symbols_list)
    molrec['real'] = [True] * len(symbols_list)
    molrec['units'] = 'Bohr'

    # Sanitize probe symbols
    molrec['symbols'] = np.asarray([symbol if symbol != 'X' else probe_symbol for symbol in molrec['symbols']])
    molrec['elem'] = molrec['symbols']

    # Generate sanitary XYZ file str 
    molstr = qcel.molparse.to_string(molrec, dtype='xyz', units='Bohr')

    return molstr

def aug2xyz(augmented: Molecule, probe_symbol: str = 'He', write: Union[str,None] = None) -> str:
    """Prepare an XYZ file defining a molecule augmented with NICS probes.

    Parameters
    ----------
    augmented : Molecule
        Molecule augmented with NICS probes for which to write to file.
    probe_symbol : {'He', 'GH', 'X'}
        Symbol to define NICS probe(s) in XYZ file. 
    write : Union[str, None], default=None
        Do the writing here?

    Returns
    ------- 
    str
        Contents of XYZ file for preservation

    Notes
    -----
    1. The appropriate probe symbol will depend on the downstream use of the XYZ file.
    Possible applications & recommended symbols:
        - Visualization with IQMol: 'He'
        - Q-Chem format: 'GH'
        - Internal/Psi4/QCElemental format: 'X'
    """
    # Sanitize molecule for writing to XYZ; see Note 1 above
    molstr = sanitize_augmented(augmented, probe_symbol=probe_symbol, dtype='xyz')

    # Build new molecule instance from temporary file.
    with tempfile.NamedTemporaryFile('w+t') as temp:
        temp.write(molstr)
        temp.seek(0)
        sanitary = qcel.models.molecule.Molecule.from_file(temp.name, dtype='xyz+')

    # Write to file
    if write is not None:
        sanitary.to_file(write, dtype='xyz')

    # Return XYZ contents to caller
    return sanitary.to_string('xyz')

def prep_g16_input(mol: Molecule, 
                   chkfile: str = 'input.chk',
                   method: str = "B3LYP",
                   basis: str = '6-311++G**',
                   title: str = '',
                   broken_symmetry: bool = False,
                   **kwargs) -> str:
    """Prepares the contents of a Gaussian16 input file for a NICS computation on molecule `mol`.

    Parameters
    ----------
    mol : Molecule
        Molecule on which to perform NICS analysis.
    chkfile : str
        Name of `.chk` file to set in G16 input
    method : str, default='B3LYP'
        Electronic structure method to use in computation
    basis : str, default='6-311++G**'
        One-electron basis set to use in computation
    broken_symmetry : {False, True}
        Break spin/spatial symmetry by mixing HOMO & LUMO?

    Returns
    -------
    nmrjob : str
        Contents of input file for computing NICS at desired level-of-theory.
    bsjob : str or None
        Contents of input file for preparing initial broken-symmetry wavefunction/density at desired
        level-of-theory.

    Raises
    ------
    Exception
        If the passed molecule `mol` does not contain NICS probes.

    Notes
    -----
    1. If a "conventional" (i.e., non-broken--symmetry) model chemistry is requested, no additional
    SCF-related options are set, so the default guess and convergence algorithm are used.
    2. If a broken-symmetry wavefunction/density is requested, the SCF options `guess=(mix,always)`
    and `scf=qc` are set
    """
    # Build comment for header
    # Check whether there are NICS probes in molecule
    if 'X' not in mol.dict()['symbols']:
        raise Exception(f"No NICS probes in molecule {mol.name}! Exiting.")

    # Build molstrings for different input files
    auglines = sanitize_augmented(mol, probe_symbol='Bq').split('\n')[2:-1]
    aug_molstr = '\n'.join(auglines)

    # Build broken symmetry job?
    if broken_symmetry:
        inpstr = fr"""%chk={chkfile}
%nproc={kwargs['nproc'] if 'nproc' in kwargs.keys() else 2}
#{method}/{basis} guess=(mix,always) scf=qc nmr=giao

{title} {mol.name}: Broken-symmetry requested; mixing HOMO/LUMO in SCF guess before NICS job

{int(mol.dict()['molecular_charge'])} {mol.dict()['molecular_multiplicity']}
{aug_molstr}

"""

    else:
        inpstr = f"""%chk={chkfile}
%nproc={kwargs['nproc'] if 'nproc' in kwargs.keys() else 2}
#{method}/{basis} nmr=giao

{title} {mol.name}: NICS job

{int(mol.dict()['molecular_charge'])} {mol.dict()['molecular_multiplicity']}
{aug_molstr}

"""

    return inpstr

def prep_qchem_input(mol: Molecule, method: str = "B3LYP", basis: str = 'cc-pVDZ', comment: str = None, **kwargs) -> str:
    """Prepares the contents of a Q-Chem input file for a NICS computation on molecule `mol`.

    Parameters
    ----------
    mol : Molecule
        Molecule on which to perform NICS analysis.
    method : str, default='B3LYP'
        Electronic structure method to use in computation
    basis : str, default='cc-pVDZ'
        One-electron basis set to use in computation
    comment : str or None, default=None
        Comment to add to the input file?

    Returns
    -------
    inpstr : str
        Contents of input file for computing NICS at desired level of theory.

    Raises
    ------
    Exception
        If the passed molecule `mol` does not contain NICS probes.
    """
    # Build comment for header
    if comment is not None:
        comment_lines = f"""$comment
{comment}
$end"""
    else:
        comment_lines = ""

    # Check whether there are NICS probes in molecule
    if 'X' not in mol.dict()['symbols']:
        raise Exception("No NICS probes in molecule! Exiting.")

    # Build molecule block
    molstr = '\n'.join(sanitize_augmented(mol, probe_symbol='GH').split('\n')[2:-1])

    chgmult = f"{int(mol.dict()['molecular_charge'])} {mol.dict()['molecular_multiplicity']}"
    molecule_lines = f"""$molecule
{chgmult}
{molstr}
$end"""
    
    # Build rem block
    ## Standard rem lines
    rem_lines = f"""$rem
JOBTYPE               NMR
METHOD                {method.upper()}
BASIS                 mixed
SCF_ALGORITHM         DIIS
PURCAR                111
SEPARATE_JK           0
LIN_K                 0
CFMM_ORDER            15
GRAIN                 1
CFMM_PRINT            2
CFMMSTAT              1
PRINT_PATH_TIME       1
LINK_MAXSHELL_NUMBER  1
SKIP_SCFMAN           0
IGUESS                core
SCF_CONVERGENCE       7
ITHRSH                10
IPRINT                23
D_SCF_CONVGUIDE       0
D_SCF_METRIC          2
D_SCF_STORAGE         50
D_SCF_RESTART         0
PRINT_PATH_TIME       1
SYM_IGNORE            1
NO_REORIENT           1
"""

    ## Additional options passed as kwargs to this function
    for option, value in kwargs:
        rem_lines += f"{option.upper()}\t\t\t{value}\n"

    ## Close rem block
    rem_lines += "$end"

    # Build basis block
    basis_lines = "$basis\n"

    for i, symbol in enumerate(mol.dict()['symbols']):
        # Recast probes as H's
        if symbol == 'X':
            symbol = 'H'
        basis_lines += f"""{symbol} {i+1}
{basis}
****
"""

    basis_lines += "$end"

    # Assemble input file
    inpstr = f"""{comment_lines}

{molecule_lines}

{rem_lines}

{basis_lines}
"""

    return inpstr

def prep_input(augmented, qcng: str, **kwargs) -> Union[str, Tuple[str]]:
    """Function to prepare input file for running NICS computation with various QC engines.

    Parameters
    ----------
    augmented : Molecule
        NICS-augmented molecule to use in input
    qcng : {'qchem', 'g16'}
        QC engine to specify input file format
    **kwargs
        Specific options to pass to the input file constructor

    Returns
    -------
    Union[str, Tuple[str]]
        String repr of the input file(s) contents, formatted for QC engine `qcng`
    
    Notes
    -----
    If G16 engine is requested and , *two* inputs will be returned
    """
    if qcng == 'qchem':
        return prep_qchem_input(augmented, **kwargs)
    
    elif qcng == 'g16':
        return prep_g16_input(augmented, **kwargs)

    else:
        raise Exception(f"Requested QC engine {qcng} not currently supported for input file generation.")

def generate_nics(xyzfile: str, 
                  ringfinder: str = 'auto',
                  xyz2molkw: Dict = {},
                  axiskw: Dict = {},
                  nicskw: Dict = {},
                  write: Union[None, List[str]] = None, 
                  inpkw: Dict = {},
                  aug2xyzkw : Dict = {},
                  wdir: Union[None, str] = None,
                  wpath: str = '.',
                  wbasename: Union[None, str] = None,
                  **kwargs
                 ) -> Dict:
    """Helper function to automate NICS generation, beginning from an unadorned XYZ file.

    Parameters
    ----------
    xyzfile : str
        Name of XYZ file to read structure from
    ringfinder : {'auto', 'manual'}
        Find rings automatically with `find_rings()`?
    xyz2molkw : Dict, default={}
        Optional kwargs to pass to `xyz2mol()`
    axiskw : Dict, default={}
        Optional kwargs to pass to `principal_inertial_axis()`
    nicskw : Dict, default={}
        Optional kwargs to pass to `generate_nics()`
    write : Union[None, List[str]], default=None
        Write generated NICS data to disk at `f"{filebase}.ext"`? Does nothing if `write=None`.
        File extension `ext` is determined by requesting one or more of the following filetypes:
            * `"xyz"` : Writes an augmented XYZ file with extension `".xyz"`
            * `"qchem"` : Writes a Q-Chem--style input file with extension `".in"`
            * `"g16"` : Writes a Gaussian 16--style input file with extension `".in"`
    inpkw : Dict, default={}
        Optional kwargs to pass to `prep_input()`
    aug2xyzkw : Dict, default={}
        Optional kwargs to pass to `aug2xyz()`
    wdir : Union[None, str], default=None
        Optional filepath to target output directory for writing
    wbasename : Union[None, str], default=None
        Manually specify base name of file to write to
    Returns
    -------
    nics : Dict
        Schema-style record of generated NICS data & provenance, with the following contents:
    """
    # Build NICS schema
    nics = {}
    nics['provenance'] = {}
    nics['provenance']['original_xyzfile'] = xyzfile

    # Read molecule structure from XYZ/XYZ+ file
    nics['provenance']['xyz2molkw'] = xyz2molkw
    nics['original_mol'] = xyz2mol(nics['provenance']['original_xyzfile'], **nics['provenance']['xyz2molkw'])

    # Automatically find rings? If not, set `rings = [None]`
    if ringfinder == 'auto':
        nics['provenance']['ringfinder'] = 'auto'
        nics['rings'] = find_rings(nics['original_mol'])
        if nics['rings'] == None:
            warnings.warn(f"WARNING: No rings detected in {nics['provenance']['original_xyzfile'].split(os.sep)[-1].split('.')[0]}. Proceed with extreme caution!!")
            nics['rings'] = [None]
    else:
        # Raise warning if user hasn't either requested ring detection or defined mol subset manually
        try:
            nics['rings'] = [axiskw['subset']]
            nics['provenance']['ringfinder'] = 'manual'
        except KeyError:
            warnings.warn("WARNING: Rings not automatically detected and no ring(s) defined manually for {nics['original_mol'].name}. Proceed with extreme caution!!")
            nics['rings'] = [None]
            nics['provenance']['ringfinder'] = None

    # TODO start here to update for using NICS schema
    # Iterate over rings, build inertial frame(s)
    nics['inertial_frames'] = []
    for ring in nics['rings']:
        system = {}
        # Solve inertial system for test molecule
        system['C'], system['axis'] = principal_inertial_axis(nics['original_mol'], subset=ring, **axiskw)
        nics['inertial_frames'].append(system)

    # Iterate over inertial frames, generate NICS for each frame
    tmpmol = nics['original_mol'].copy() 
    for frame in nics['inertial_frames']:
        # Generate NICS from centroid along axis
        tmpmol = augment(tmpmol, frame['C'], frame['axis'], **nicskw)
    nics['augmented'] = tmpmol

    # Write data?
    if write is None:
        nics['files'] = None
        return nics

    else:
        # Check for manually provided basename
        if wbasename is not None:
            filebase = wbasename
        else:
            # Get filebase for input file writing
            filebase = nics['provenance']['original_xyzfile'].split(os.sep)[-1].split('.')[0]

        if isinstance(write, str):
            write = [write]

        # Assemble files to write
        nics['files'] = {}
        for wtype in write:
            if wtype in ['qchem', 'g16']:
                # Build MC?
                try:
                    method = inpkw['method']
                    basis = inpkw['basis'].replace('*','s').replace('+','p')
                    fname = f"{filebase}_geom-{'BS' if inpkw['broken_symmetry'] else ''}{method}_{basis}_nics"
                except KeyError:
                    warnings.warn("No method or basis set name passed. Generic filename it is!")
                    fname = f"{filebase}_geom-nics"
            if wtype == 'qchem':
                # Generate input content once all the probes have been placed
                contents = prep_qchem_input(nics['augmented'], **inpkw)

                nics['files'][wtype] = {}
                nics['files'][wtype]['dirpath'] = wpath
                nics['files'][wtype]['name'] = f"{fname}.in"
                nics['files'][wtype]['contents'] = contents

            if wtype == 'g16':
                # Make sure broken symmetry kwarg set
                try:
                    testing = inpkw['broken_symmetry']
                except KeyError:
                    inpkw['broken_symmetry'] = False

                # Make sure checkfile & extension is set properly if it's g16
                inpkw['chkfile'] = f"{fname}.chk"

                # Generate input content once all the probes have been placed
                contents = prep_g16_input(nics['augmented'], **inpkw)

                nics['files'][wtype] = {}

                nics['files'][wtype] = {}
                nics['files'][wtype]['dirpath'] = wpath
                nics['files'][wtype]['name'] = f"{fname}.com"
                nics['files'][wtype]['contents'] = contents

            elif wtype == 'xyz':
                # Generate XYZ contents
                contents = aug2xyz(nics['augmented'], **aug2xyzkw)

                nics['files']['xyz'] = {}
                nics['files']['xyz']['name'] = f"{filebase}-nicsviz.xyz"
                nics['files']['xyz']['dirpath'] = wpath
                nics['files']['xyz']['contents'] = contents

            else:
                raise Exception(f"Requested write format {wtype} not supported.")

        # Write the files
        for wtype, details in nics['files'].items():
            fullname = details['dirpath'] + os.sep + details['name']
            contents = details['contents']
            with open(fullname, 'w+') as f:
                f.write(contents)

        return nics
    
def which_engine(output: str) -> str:
    """Determines which qcng produced an output file containing NICS data.

    Parameters
    ----------
    output : str
        Name of output file whose engine to determine

    Returns
    -------
    qcng : str
        Name of qcng which produced `output`.

    Raises
    ------
    Exception
        If the generating QC engine for the output file `output` cannot be determined

    Notes
    -----
    1. Currently, only Gaussian and Q-Chem engines are recognized; other output
    file types will raise an exception.
    """

    with open(output, 'r') as f:
        lines = f.readlines()

    # Iterate over lines, return if keyword found
    for line in lines[:100]:
        if 'Gaussian' in line: return 'g16'

        if 'Q-Chem' in line: return 'qchem'

        if 'Psi4' in line: return 'psi4'

    raise Exception(f"Generating QC engine for output file {output} could not be determined! Are you sure it's an output file?")

def extract_nics(output: str, write: Union[None, str] = None) -> Dict:
    """Extracts NICS data from output file `output`.

    Parameters
    ----------
    output : str
        Name of output file from which to extract NICS data
    write : Union[None, str], default=None
        Write extracted NICS data to augmented XYZ file?

    Returns
    -------
    nics : Dict
        Dict containing record of NICS computation & NICS data

    Raises
    ------
    Exception
        If output file is determined to be an unsupported type
    """

    # Which QC engine generated the output?
    qcng = which_engine(output)

    # Extract NICS data from output
    if qcng == 'g16':
        nics = parse_g16_output(output)

    elif qcng == 'qchem':
        nics = parse_qchem_output(output)

    else:
        raise Exception(f"Output file {output} was generated by QC engine {qcng}, which is currently not supported! Exiting.")

    # Store extracted data to file?
    if write is not None:
        nics2xyz(nics, write)

    return nics

def parse_g16_output(output: str) -> Dict:
    """Extracts NICS data and molecule from Gaussian16 output.

    Parameters
    ----------
    output : str
        Output file (typically with extension `.log`) from G16 NICS job

    Returns
    -------
    nics : Dict
        Dict repr of molecule & NICS analysis
    """
    # Dict setup
    nics = {}
    nics['job'] = {}
    nics['molecule'] = {}

    # Get name of file for posterity
    nics['job']['output'] = output

    # First things first: read the file
    with open(output, 'r') as f:
        outlines = f.readlines()

    # Iterate over lines, save important stuff as we go
    molstart = molend = None
    isotropic = {}
    for line in outlines:
        if 'nproc' in line:
            nics['job']['nproc'] = int(line.split('=')[-1].strip())
        if '#' in line and r'\\' not in line:
            nics['job']['route'] = line
        if "Input orientation:" in line:
            molstart = outlines.index(line) + 5
        if "Distance matrix" in line:
            molend = outlines.index(line) - 1
        if "Anisotropy" in line:
            stuff = line.split()
            isotropic[int(stuff[0])] = -float(stuff[4])

    # Go back and pull ranges of lines
    mollines = outlines[molstart:molend]

    moldict = {}
    for atom in mollines:
        # Handle molecular structure & probes
        stuff = atom.split()
        i = int(stuff[0])
        moldict[i] = {}
        moldict[i]['E'] = qcel.periodictable.to_E(int(stuff[1]))
        moldict[i]['x'] = float(stuff[3])
        moldict[i]['y'] = float(stuff[4])
        moldict[i]['z'] = float(stuff[5])

        # Get isotropic deshielding value, add to moldict
        moldict[i]['-iso'] = isotropic[i]

    nics['molecule'] = moldict

    return nics

def parse_qchem_output(output: str) -> Dict:
    raise Exception("Q-Chem NICS output parsing not yet implemented!")

def nics2xyz(nics: Dict, inplace: bool = False, name: Union[str, None] = None) -> str:
    """Prepare extracted NICS data as augmented XYZ file, and write if indicated.

    Parameters
    ----------
    nics : Dict
        Dictionary containing extracted NICS data to translate to XYZ format
    inplace : {False, True}
        Do the writing directly here, or only return XYZ file contents?
    name : Union[str, None], default=None
        Name to manually write file to

    Returns
    -------
    xyzlines : str
        XYZ file contents containing extracted NICS data

    Notes
    -----
    1. If no xyz filename is specified manually by setting the `name` kwarg, the filename
    will be inferred from the name of the output file from which the data were originally
    extracted. Furthermore, the XYZ file will be written alongside that output file in
    the directory tree.
    """
    # Begin XYZ file
    xyzlines = f"""{len(nics['molecule'])}
Data extracted from output file: {nics['job']['output']}
"""
    # write atoms
    for atom in nics['molecule'].values():
        if atom['E'] != 'X':
            xyzlines += f"{atom['E']:<4s} {atom['x']:>13.10f} {atom['y']:>13.10f} {atom['z']:>13.10f}\n"
        else:
            xyzlines += f"{atom['E']:<4s} {atom['x']:>13.10f} {atom['y']:>13.10f} {atom['z']:>13.10f} {atom['-iso']:>13.10f}\n"
            
    # Do writing here?
    if inplace:
        # Get name from NICS output file? Or is it set manually?
        if name is not None:
            fname = name
        else:
            fname = os.path.splitext(nics['job']['output'])[0] + '.xyz'

        with open(fname, 'w+') as f:
            f.write(xyzlines)

    return xyzlines

# ==> Testing <==

def test_monocyclic():
    # Read test molecule
    mol = xyz2mol('test/3cTL.xyz', molecular_charge=2, molecular_multiplicity=1)

    # Solve inertial system for test molecule
    C, axis = principal_inertial_axis(mol, subset=[0,1,2,3,4,5,6,7], masswt=False)

    # Generate NICS from centroid along axis, save to xyz
    augmented = augment(mol, C, axis, recipe='IsoSimple')

    inptest = prep_qchem_input(augmented)

    # Write input to {filebase}.in
    with open(f"test/3cTL-nics.in", 'w+') as f:
        f.write(inptest)
    
def test_polycyclic():
    # Build 3d from xyx
    mol = xyz2mol('test/3d.xyz', molecular_charge=2, molecular_multiplicity=1)

    # Find rings in 3d
    rings = find_rings(mol)
    
    # Build inertial systems for both rings
    systems = []
    for ring in rings:
        system = {}
        C, axis = principal_inertial_axis(mol, subset=ring, masswt=False)
        system['C'] = C
        system['axis'] = axis
        systems.append(system)

    # Augment mol with probes for each system
    augmented = mol.copy()
    for system in systems:
        augmented = augment(augmented, system['C'], system['axis'], recipe='IsoSimple') 

    # Write aug to XYZ for viz & verification
    _ = aug2xyz(augmented, 'test/3d-nics.xyz', probe_symbol='He')

def test_2c():
    # Build 3d from xyx
    mol = xyz2mol('test/2c.xyz', molecular_charge=1, molecular_multiplicity=1)

    # Find rings in 3d
    rings = find_rings(mol)
    
    # Build inertial systems for both rings
    systems = []
    for ring in rings:
        system = {}
        C, axis = principal_inertial_axis(mol, subset=ring, masswt=False)
        system['C'] = C
        system['axis'] = axis
        systems.append(system)

    # Augment mol with probes for each system
    augmented = mol.copy()
    for system in systems:
        augmented = augment(augmented, system['C'], system['axis'], recipe='IsoSimple') 

    # Write aug to XYZ for viz & verification
    _ = aug2xyz(augmented, write='test/viz/2c-nics.xyz', probe_symbol='He')

def test_driver():
    # Test both mono & polycyclic
    #geoms = ['test/3cTL.xyz', 'test/3d.xyz']
    geoms = ['test/3cTL.xyz']

    for xyz in geoms:
        nicskw = {'recipe': 'IsoSimple'}
        nics = nics_driver(xyz, nicskw=nicskw)
        #print(nics)

    return True

def test_g16_input():
    # Read test molecule
    mol = xyz2mol('test/2c.xyz', molecular_charge=1, molecular_multiplicity=1)
    fbase = '2c'

    # Solve inertial system for test molecule
    C, axis = principal_inertial_axis(mol, subset=[0,1,2,3,4,5,6,7], masswt=False)

    # Generate NICS from centroid along axis, save to xyz
    augmented = augment(mol, C, axis, recipe='IsoSimple')

    inptest = prep_g16_input(augmented, method="UB3LYP", chkfile=f"{fbase}.chk", title='nics test')

    return inptest
    

# ==> Main <==
if __name__ == "__main__":
    # Parser setup
    parser = argparse.ArgumentParser(description="Augment molecular geometries with NICS probes & prepare input files for NICS computation.")
    # Testing
    parser.add_argument('-t', '--test', action="store_true", help="Quick test of NICS helpers")
    # XYZ -> Molecule
    parser.add_argument('-g', '--geometry', default=None, help="Name of XYZ file(s) to augment with NICS probes", nargs='*')
    parser.add_argument('-c', '--chgmult', default=[0, 1], help="Molecular charge and multiplicity for geometry", nargs=2)
    # NICS probe recipes
    parser.add_argument('-r', '--recipe', default='IsoSimple', help="NICS recipe(s) to follow", nargs='*')
    # Write data?
    parser.add_argument('-f', '--filetype', default=None, choices=['xyz', 'qchem', 'g16'], help="File types to write augmented molecules to", nargs='*')
    # Input file options
    parser.add_argument('-m', '--modelchem', default='B3LYP/cc-pVDZ', help="Combination of method/basis to use when preparing input file(s)")
    # Write locations
    parser.add_argument('-l', '--location', default='.', help="Path to directory within which to write files")
    # Optional filename
    parser.add_argument('-n', '--name', default=None, help="Manually specify output filename")
    # Parse NICS output?
    parser.add_argument('-e', '--extract', default=None, help="Output file(s) from which to extract NICS values", nargs='*')
    
    args = parser.parse_args()

    # Testing: manually find 3d rings & prep augmented mol
    if args.test:
        assert test_driver(), "Shit's fucked!"

        # Testing g16
        #print(test_g16_input())
        sys.exit()

    # Running to parse output(s)?
    if args.extract is not None:
        # Handle files to extract
        if isinstance(args.extract, list):
            outs = args.extract
        else:
            outs = glob.glob(args.extract)
    
        # Iterate over output files
        for out in outs:
            # Parse output, return NICS data
            nics = extract_nics(out)

            # Write NICS data to XYZ file named in same way as output
            nics2xyz(nics, inplace=True)

        # If we're running to extract data, that's all we're doing. 
        sys.exit()

    # Listify options if they're a single value
    # Collect specified geometries
    if isinstance(args.geometry, list):
        geoms = args.geometry
    else:
        geoms = glob.glob(args.geometry)
    
    # Handle recipes
    if isinstance(args.recipe, list):
        recipes = args.recipe
    else:
        recipes = [args.recipe]

    # Handle writes
    if isinstance(args.filetype, list):
        outs = args.filetype
    else:
        outs = [args.filetype]

    # Handle modelchem
    method = args.modelchem.split('/')[0]
    basis = args.modelchem.split('/')[1]

    # Wanting to do broken symmetry?
    if 'BS-' in method:
        # Remove prefix from method
        method = method.replace('BS-', '')
        # Toggle broken_symmetry
        broken_symmetry=True
    else:
        broken_symmetry=False
    
    # Now do work!
    for geom in geoms:
        # Iterate over NICS recipes
        for recipe in recipes:
            # Iterate over requested write types
            for write in outs:
                nics = generate_nics(geom,
                            ringfinder='auto',
                            xyz2molkw={'molecular_charge':int(args.chgmult[0]), 'molecular_multiplicity':int(args.chgmult[1])},
                            axiskw={},
                            nicskw={'recipe': recipe},
                            inpkw={'method': method, 'basis': basis, 'broken_symmetry': broken_symmetry},
                            write=args.filetype,
                            wpath=args.location,
                            wbasename=args.name
                        )



