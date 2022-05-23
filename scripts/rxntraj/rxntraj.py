"""Builds a straight-line reaction trajectory via structural interpolation."""

__author__ = "Dominic A. Sirianni, PhD"
__license__ = "BSD-3-Clause"
__copyright__ = "(c) 2020-2022 The Parish Lab"
__date__ = "2021-01-13"

import os
import psi4
import numpy as np
import scipy as sp
from copy import deepcopy
import argparse

from typing import List, Dict, Tuple, Optional, Iterable, Callable, Any, Union

def xyz2mol(xyzfile: str) -> Dict[str, Any]:
    """Reads geometry from XYZ file, returns dict w/ mol info

    """
    mol = {}
    # Use NumPy loadtxt to directly load molecule
    mol['atoms'] = np.loadtxt(xyzfile, skiprows=2, usecols=0, dtype=str)
    mol['coords'] = np.loadtxt(xyzfile, skiprows=2, usecols=(1,2,3))

    return mol

def xyz2p4mol(xyzfile: str, units: str ='AA') -> Any:
    """Reads geometry from XYZ file, returns psi4.core.Molecule

    """
    # Get geometry from string of concatenated XYZ lines
    with open(xyzfile, 'r') as xyz: 
        lines = xyz.readlines()[2:]
    mol = psi4.geometry(''.join(lines))
    if units == 'AA':
        mol.set_units(psi4.core.GeometryUnits.Angstrom)
    elif units == 'a0':
        # Default for QChem
        mol.set_units(psi4.core.GeometryUnits.Bohr)
    else:
        raise Exception(f"XYZ coordinates must be declared in either AA or a0, not {units}.")

    return mol

def get_normal_modes(mol: Any) -> Dict[Union[float, complex], Any]:
    """Gets normal modes of vibration for molecule `mol` using Psi4.

    """
    # Send psi4 output to file
    psi4.core.set_output_file('output.dat', False)

    # Compute molecular Hessian, get basic quantities
    e, wfn = psi4.frequencies('hf/sto-3g', molecule=mol, return_wfn=True)
    natoms = mol.natom()
    masses = np.array([mol.mass(i) * psi4.constants.conversion_factor("amu", "atomic_unit_of_mass") for i in range(mol.natom())])

    # Mass-weight Hessian
    H = np.asarray(wfn.hessian())
    M = np.diag(1 / np.sqrt(np.repeat(masses, 3))) # Matrix of masses for weighting
    mH = M.T.dot(H).dot(M)

    # Diagonalize Hessian
    k2, Lxm = np.linalg.eigh(mH)

    # Get freqs from ew of Hessian
    freqs = []
    mode = 3 * natoms - 1 # Last vibrational mode
    while mode >= 6: # Only keep vib modes, not rotation or translation
        # Real vib mode?
        if k2[mode] >= 0.0:
            freqs.append(np.sqrt(k2[mode]) * psi4.constants.conversion_factor("hartree^(1/2)/(bohr * atomic_unit_of_mass^(1/2))", "cm^-1") / (2 * np.pi))
        # Imaginary vib modes get a j
        else:
            freqs.append((np.sqrt(abs(k2[mode])) * psi4.constants.conversion_factor("hartree^(1/2)/(bohr * atomic_unit_of_mass^(1/2))", "cm^-1") / (2 * np.pi))*j)
        mode -= 1

    # Un mass-weight mode vectors (m_e^{-1/2} -> unitless)
    Lx = M.dot(Lxm)
    # Reshape mode vectors (1,natoms) -> (natoms, 3), zip with vib freqs
    modes = {}
    for i, freq in enumerate(freqs):
        modes[freq] = Lx[:,6+i].reshape(natoms,3)

    return modes
    
def get_rxn_modes(start: Any, end: Any, thresh=None) -> Tuple[Any, Dict[str, Any]]:
    """Writes rxn displacement vector `xi` in basis of vib modes of molecule `start`.

    Parameters
    ----------
    start : psi4.core.Molecule
        Reactant molecule
    end : psi4.core.Molecule
        Product molecule
    thresh : float64, optional
        Threshold to prune low-importance vib modes from xi

    Returns
    -------
    xi : Dict[float, Any]
        {component: mode} dict containing decomposition of `xi` in normal mode basis. 
    """
    # Build xi vector
    xi = rxnpath(start, end, usepsi=True, onlyxi=True)
    
    # Get normal modes of `start`
    modes = get_normal_modes(start)

    # Get components of xi in basis of normal modes
    components = np.zeros(len(modes.keys()))
    for nu, mode in modes.items():
        # Get index for mode (modes in decreasing order)
        i = list(modes.keys()).index(nu)
        # < mode | xi > => np.einsum('ij,ij->', mode, xi)
        components[i] = np.einsum('ij,ij->', mode, xi)

    # Return xi vector components & vib modes
    return components, modes

def rxnpath(start: Any, end: Any, intervals: int = 1, usepsi=False, onlyxi=False) -> Dict[float, Any]:
    """Constructs a reaction path between molecule `start` and `end`,
    with `intervals` many interpolated molecules along reaction pathway xi.

    """
    if usepsi:
        # Check that the atoms are in the same bloody order. You're fucked if they're not.
        if False in (start.to_arrays()[-1] == end.to_arrays()[-1]):
            raise Exception("Bloody atoms aren't in the same order in start/end.")
        # Align molecules
        align = end.B787(start, run_resorting=True, run_to_completion=True, run_mirror=True)
        end = align[-1]
        
        # Send mols to arrays, build traj from there
        # TODO incorporate alignment
        startarr = start.to_arrays()
        endarr = end.to_arrays()

        # Build displacement vector from start -> end
        xi = endarr[0] - startarr[0]
        
        if onlyxi:
            return xi

        traj = {}
        for disp in np.arange(0, 1 + 1/intervals, 1/intervals):
            # Step `disp` far along xi
            newcoords = startarr[0] + disp * xi
            newcoords = psi4.core.Matrix.from_array(newcoords)
            
            # Cast newcoords -> psi4.core.Matrix, use to create new mol w/ new coords
            newmol = start.clone()
            newmol.set_geometry(newcoords)
            newmol.update_geometry()
            traj[disp] = newmol
    else:
        # Check that the atoms are in the same bloody order. You're fucked if they're not.
        if False in (start['atoms'] == end['atoms']):
            raise Exception("Bloody atoms aren't in the same order in start/end.")
        if np.allclose(np.array(start['coords']), np.array(end['coords'])):
            raise Exception("Geometry is the same!!")
        # No alignment, you better know what you're doing!
        startarr = start['coords']
        endarr = end['coords']

        # Build displacement vector from start -> end
        xi = endarr - startarr

        if onlyxi:
            return xi

        traj = {}
        for disp in np.arange(0, 1 + 1/intervals, 1/intervals):
            # Step `disp` far along xi
            newcoords = startarr + disp * xi
            
            # Cast newcoords -> psi4.core.Matrix, use to create new mol w/ new coords
            newmol = deepcopy(start)
            newmol['coords'] = newcoords
            traj[disp] = newmol

    return traj

def save_rxn_traj(path: Dict[float, Any], filename: str = 'traj', filepath: str = '.', gen_input=False) -> None:
    """Saves all displaced geometries along reaction path in XYZ format.

    """
    for disp, mol in path.items():
        d = np.around(np.array(disp), decimals=2)
        if gen_input:
            pass
        else:
            if isinstance(mol, dict):
                with open(f"{filepath}{os.sep}{filename}_{d}.xyz", 'w+') as f:
                    f.write(f"{len(mol['atoms'])}\n\n")
                    for atom, coords in zip(mol['atoms'], mol['coords']):
                        f.write(f"{atom:<2} {coords[0]:>10.2f} {coords[1]:>10.2f} {coords[2]:>10.2f}\n")
            else:
                mol.save_xyz_file(f"{filepath}{os.sep}{filename}_{d}.xyz", True)

# TODO: Input writer for snapshots

def prep_qchem_input(mol: Molecule,
                     template: str = None,
                     method: str = None,
                     basis: str = None,
                     comment: str = None,
                     **kwargs
                    ) -> str:
    """Prepares the contents of a Q-Chem input file for a computation on molecule `mol`.

    Parameters
    ----------
    mol : Molecule
        Molecule on which to perform NICS analysis.
    template : str or None, default=None
        Template input file contents to replicate with molecule `mol`

    Returns
    -------
    inpstr : str
        Contents of input file.
    """
    # Working from a template?
    if template is not None:
        # Read template, join lines
        with open(template, 'r') as f:
            lines = template.readlines()

        inpstr = '\n'.join(lines)

    # Build molecule block
    molstr = '\n'.join(mol.to_string(dtype='xyz').split('\n')[2:-1])

    chgmult = f"{int(mol.dict()['molecular_charge'])} {mol.dict()['molecular_multiplicity']}"
    molecule_lines = f"""$molecule
{chgmult}
{molstr}
$end"""

    return inpstr.format(chgmult,molstr)

def prep_input(mol: Molecule, template: str, **kwargs) -> str:
    """Function to prepare input file for running NICS computation with various QC engines.

    Parameters
    ----------
    mol : Molecule
        Molecule to use in generated input
    template : str
        Name of file from which to pull input file template
    **kwargs
        Specific options to pass to the input file constructor

    Returns
    -------
    str
        String repr of the input file(s) contents, formatted for QC engine `qcng`
    """
    # Determine qcng
    qcng = which_engine(template)

    if qcng == 'qchem':
        contents = parse_qchem_template(template)
        return prep_qchem_input(mol, template=contents, **kwargs)
    
    elif qcng == 'g16':
        return prep_g16_input(mol, **kwargs)

    elif qcng == 'psi4':
        return prep_psi4_input(mol, **kwargs)

    else:
        raise Exception(f"Requested QC engine {qcng} not currently supported for input file generation.")

def which_engine(inpfile: str) -> str:
    """Determines which qcng a given template input file is written for

    Parameters
    ----------
    inpfile : str
        Name of output file whose engine to determine

    Returns
    -------
    qcng : str
        Name of qcng for which `inpfile` is formatted.

    Raises
    ------
    Exception
        If the QC engine corresponding to the `inpfile` format cannot be determined

    Notes
    -----
    1. Currently, only Q-Chem, Psi4, and Gaussian-formatted template input files are recognized.
    """

    with open(output, 'r') as f:
        lines = f.readlines()

    # Iterate over lines, return if keyword found
    for line in lines[:100]:
        if '%chk' in line: return 'g16'

        if '$molecule' in line: return 'qchem'

        if 'molecule' and '{' in line.split(): return 'psi4'

    raise Exception(f"QC engine for template input file {inpfile} could not be determined! Are you sure it's formatted correctly?")

# ==> Testing <==
if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Creates a reaction path between two molecular endpoints.")
    parser.add_argument('start', help="Name of XYZ file containing starting molecular geometry")
    parser.add_argument('end', help="Name of XYZ file containing ending molecular geometry")
    parser.add_argument('increment', help="Number of increments to create along reaction path")
    parser.add_argument('-p', '--path', default='.', help="File path to use for XYZ structure files along reaction trajectory.")
    parser.add_argument('-n', '--name', default='traj', help="File name to use for all XYZ structure files along reaction trajectory.")
    parser.add_argument('-e', '--engine', default='np', help="Engine to use for mol handling.")
    parser.add_argument('-t', '--test', default=False, help="Toggles testing")
    args = parser.parse_args()

    if args.test:
        pass
    else:
        if args.engine == 'psi':
            start = xyz2p4mol(args.start)
            end = xyz2p4mol(args.end)
            traj = rxnpath(start, end, int(args.increment), usepsi=True)
            save_rxn_traj(traj, filename=args.name, filepath=args.path)
    
        else:
            start = xyz2mol(args.start)
            end = xyz2mol(args.end)
            traj = rxnpath(start, end, int(args.increment))
            save_rxn_traj(traj, filename=args.name, filepath=args.path)
