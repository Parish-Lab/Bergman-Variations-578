"""
iwPert is a short script for extracting vibrational modes from a Q-Chem
frequency analysis, and perturbing a geometry along one or more desired
vibrational modes.
"""

__authors__ = "Dominic A. Sirianni"
__license__ = "BSD-3-Clause"
__copyright__ = "(c) 2020-2022, The Parish Lab"
__date__ = "2020-11-19"

import numpy as np
from io import StringIO
import collections.abc
import itertools as it

import glob
import argparse
from typing import List, Dict, Tuple, Optional, Iterable, Callable, Any, Union

def load_output(filepath: str) -> List[str]:
    """Loads output file at `filepath` into memory.

    Parameters
    ----------
    filepath : str
        Path from cwd to desired output file for scraping

    Returns
    -------
    outlines : List[str]
        List of lines in the output file `filepath`
    """

    with open(filepath, 'r') as f:
        outlines = f.readlines()

    return outlines

def extract_vib_analysis(outlines: List[str]) -> Dict[int, Dict[str, Any]]:
    """Extracts vibrational analysis from Q-Chem output.

    Parameters
    ----------
    outlines : List[str]
        Lines from output file, loaded as a list from `load_output()`

    Returns
    -------
    modes : Dict[int: Dict[str: Any]]
        Dict of dicts corresponding to vibrational mode info.
    """
    # Get start & end info, extract vib lines
    for l in outlines:
        # Find start of vib analysis
        if "VIBRATIONAL ANALYSIS" in l:
            start = outlines.index(l) + 10
        # Find end of vib analysis
        elif "STANDARD THERMODYNAMIC QUANTITIES AT" in l:
            end = outlines.index(l)
    viblines = outlines[start:end]

    # Separate viblines into paragraph blocks
    blocks = ''.join(viblines).split('\n\n')[:-1]
    
    modes = {}
    # Iterate over blocks, split into modes w/in blocks
    for b in blocks:
        # Extract block info
        blines = b.split('\n')
        bmodes = [int(i) for i in blines[0].split()[1:]]
        freqs = [np.inf if '**' in i else float(i) for i in blines[1].split()[1:]]
        ks = [np.inf if '**' in i else float(i) for i in blines[2].split()[2:]]
        us = [np.inf if '**' in i else float(i) for i in blines[3].split()[2:]]
        IRactive = [True if i == 'YES' else False for i in blines[4].split()[2:]]
        IRintens = [np.inf if '**' in i else float(i) for i in blines[5].split()[2:]]
        Raman = [True if i == 'YES' else False for i in blines[6].split()[2:]]
        nmodes = len(bmodes)

        # Treat atom:mode lines like numpy array to split
        linesforatoms = StringIO('\n'.join(blines[8:-1]))
        atoms = np.loadtxt(linesforatoms, dtype=str, usecols=0)
        disps = []
        disps.append(np.loadtxt(StringIO('\n'.join(blines[8:-1])), usecols=(1,2,3)))
        disps.append(np.loadtxt(StringIO('\n'.join(blines[8:-1])), usecols=(4,5,6)))
        disps.append(np.loadtxt(StringIO('\n'.join(blines[8:-1])), usecols=(7,8,9)))
        
        # Do same treatment for TransDip        
        if len(blines[-1]) != 10:
            transdisps = [None, None, None]
        else:
            transdisps = []
            transdisps.append(np.loadtxt(StringIO(blines[-1]), usecols=(1,2,3)))
            transdisps.append(np.loadtxt(StringIO(blines[-1]), usecols=(4,5,6)))
            transdisps.append(np.loadtxt(StringIO(blines[-1]), usecols=(7,8,9)))

        # Store block info 
        for m in range(nmodes):
            mode = {}
            mode['Freq'] = freqs[m]
            mode['Force Const'] = ks[m]
            mode['Red. Mass'] = us[m]
            mode['IR Active'] = IRactive[m]
            mode['IR Intens'] = IRintens[m]
            mode['Raman Active'] = Raman[m]
            mode['atoms'] = atoms
            mode['disp'] = disps[m]
            mode['TransDip'] = transdisps[m]
            modes[bmodes[m]] = mode

    return modes
    
def extract_geom(outlines: List[str]) -> Dict[str, Any]:
    """Extract final molecular geometry from output in `outlines`.

    Parameters
    ----------
    outlines : List[str]
        Lines of output file to parse for molecule

    Returns
    -------
    molecule : Dict[str, Any]
        Dict w/ atoms & XYZ coordinates for molecule
    """
    # Pull final geom lines from output
    geom = []
    for l in outlines[::-1]:
        if "Nuclear Orientation" in l:
            start = outlines.index(l) + 3
    for l in outlines[start:]:
        if '--------------------' in l:
            break
        else:
            geom.append(l)
    
    # Read geom to arrays with numpy
    molstr = '\n'.join(geom)
    atoms = np.loadtxt(StringIO(molstr), dtype=str, usecols=1)
    xyz = np.loadtxt(StringIO(molstr), usecols=(2, 3, 4))

    # Package atoms & XYZ into molecule dict
    molecule = {'atoms': atoms, 'coords': xyz}

    return molecule

def perturb(mol: Dict[str, Any], \
            modes: Union[Dict[str, Any], Dict[float, Dict[str, Any]]], \
            mag: Optional[Union[float, int]] = 0.5) -> \
            Dict[str, Any]:
    """Perturbs molecule `mol` along vib mode(s) `modes` with magnitude `mag`.

    Perturbation of a molecular geometry G along vibrational mode V with 
    scalar magnitude n to yield a new geometry M is done according to the
    expression
        M = G + \sum_i n * V_i,
    where * represents elementwise scalar multiplication and + represents
    elementwise array addition.

    Parameters
    ----------
    mol : Dict[str, Any]
        Dict defining molecule's atoms & geometry
    modes : Union[Dict[str, Any], Dict[float, Dict[str, Any]]]
        Dict of mode along which to perturb `mol`
    mag : Optional[Union[float, int]]
        Magnitude(s) by which to perturb the original geometry. 

    Returns
    -------
    Dict[str, Any]
        Molecule dictionary containing perturbed geometry
    """
    # Check number of modes to perturb along, adjust accordingly
    if not isinstance(modes, collections.abc.Iterable):
        modes = [modes]

    # Iterate over modes, accumulate into new geometry
    newmol = mol.copy()
    for mode in modes:
        # Check to see the atomic ordering is the same for mol & mode
        molatoms = list(mol['atoms'])
        modeatoms = list(mode['atoms'])
        assert molatoms==modeatoms, f"Atomic ordering in molecule and vibrational mode do not match. Mol: {mol['atoms']}, Mode: {mode['atoms']}"

        # Check to see what caller wants re: number of perturbations
        newmol['coords'] += mag * mode['disp']

    return mol

def moldict2xyz(mol: Dict[str, Any], xyzname: Optional[str] = None, \
               comment: Optional[str] = None) -> str:
    """Write a molecule dictionary to XYZ file called `xyzname`

    Parameters
    ----------
    mol : Dict[str, Any]
        Molecule dictionary containing atomic symbols and Cartesian coordinates
    xyzname : str, Optional
        File name (and relative path location, if appropriate) to write XYZ file
        Default: None
    comment : str, Optional
        Comment to be included on the second line of the XYZ file

    Returns
    -------
    xyzlines : str
        Lines of XYZ file        
    """
    xyzlines = """"""

    # Get # of atoms
    xyzlines += f"{len(mol['atoms'])}\n"
    
    # Add comment if there is one
    if comment is not None:
        xyzlines += f"{comment}\n"
    else:
        xyzlines += "\n"
    
    # Iterate over atoms & coords, zip, append
    for i, j in zip(mol['atoms'], mol['coords']):
        atomline = f"{i}\t"
        for k in j:
            atomline += f"{k:>16.12f}\t"
        xyzlines += atomline + "\n"

    # Save XYZ?
    if xyzname is not None:
        with open(xyzname, 'w+') as xyz:
            xyz.write(xyzlines)
    
    return xyzlines

# ==> Main Execution <==
if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Extract & perturb molecules along desired vibrational modes from output of Q-Chem frequency analysis.")
    parser.add_argument('outfile', help="Name of output file to parse.")
    parser.add_argument('-m', '--mode', help="Vibrational mode along which to perturb geometry", nargs='*')
    parser.add_argument('-d', '--disp', help="Scalar magnitude to displace along mode", nargs='*')
    #parser.add_argument('-d', '--disp', help="Scalar magnitude to displace along mode")
    parser.add_argument('-w', '--writexyz', default=None, help="Location to write XYZ file(s) of perturbed geometry(ies). Filename will be appended with mode and displacement information; this argument should be considered a filename prefix.")
    args = parser.parse_args()

    ## Testing
    # Load outfile, extract vib analysis & mol
    ol = load_output(args.outfile)
    vib = extract_vib_analysis(ol)
    mol = extract_geom(ol)

    # Perform multiple displacements?
    if not isinstance(args.disp, list):
        disps = [float(args.disp)]
    else:
        disps = [float(i) for i in args.disp]

    # Iterate over displacements
    for disp in disps:
        # Cast provided mode(s) to list for simultaneous projection
        modes = []
        if not isinstance(args.mode, list):
            mode = [int(args.mode)]
        else:
            mode = [int(i) for i in args.mode]

        description = ""
        for m in mode:
            modes.append(vib[m])
            description += f"+ w_{m} * {disp} "

        # Perturb mol along indicated mode(s)
        newmol = perturb(mol, modes, mag=float(disp))

        # Build comment line for XYZ file indicating what was done
        comment = f"Perturbed geometry P = G_0 {description}"
        xyz = moldict2xyz(mol, comment=comment)

        # Write perturbed XYZ to file?
        if args.writexyz is not None:
            filename = f"{args.writexyz}_pert_{''.join(description.split())}.xyz"
            with open(filename, 'w+') as f:
                f.write(xyz)
        else:
            print(xyz)

