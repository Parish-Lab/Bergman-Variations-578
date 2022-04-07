Variations on the Bergman Cyclization Theme: Electrocyclizations of Penta-, Hepta-, and Octa-Diynes
===================================================================================================

Welcome to the supplementary repository for the 5-7-8 project! Here, we make
available the software that was written to facilitate the automatic generation,
performance, and analysis of nucleus-independent chemical shift computations
during the course of this project. 

Enjoy!

~ The Authors


## Using the `nics-helper.py` Script

### Installation

To use the helper script `nics-helper.py`, first clone this repository to
your local machine before proceeding to install the required dependencies.

It is recommended that you perform this installation in a new Conda environment
to protect against version clashes with other Python packages you may have installed
on your machine. To do this, first download and install [Anaconda]() or [Miniconda](),
before executing the following from the command line:

```bash
user@host:~$ conda create -n nicshelper python=3.9 qcelemental numpy argparse tempfile networkx warnings inspect glob copy -c conda-forge
```

This command will create a new Conda environment named `nicshelper`, into which
all the dependencies will be installed. To then use the script, first activate the
Conda environment you just created:

```bash
user@host:~$ conda activate nicshelper
```

And voila! Now you're ready to use the script.

### Examples

The `nics_helper.py` script can be used either by importing it as a Python
module into another script or by calling it directly from the command line.
Below you'll find some typical use cases and how to call the script from the
command line to perform those actions (using representative sample files 
available in the `samples/` directory).

#### Example 00: Getting a help summary from the script

To get a summary of command-line arguments to `nics-helper.py`, call the script with
the argument `-h` or `--help`:

```bash
(nicshelper) user@host:~/path/to/repo$ python nics-helper.py --help

usage: nics-helper.py [-h] [-t] [-g [GEOMETRY ...]] [-c CHGMULT CHGMULT] [-r [RECIPE ...]]
                      [-f [{xyz,qchem,g16} ...]] [-m MODELCHEM] [-l LOCATION] [-n NAME]
                      [-e [EXTRACT ...]]

Augment molecular geometries with NICS probes & prepare input files for NICS computation.

optional arguments:
  -h, --help            show this help message and exit
  -t, --test            Quick test of NICS helpers
  -g [GEOMETRY ...], --geometry [GEOMETRY ...]
                        Name of XYZ file(s) to augment with NICS probes
  -c CHGMULT CHGMULT, --chgmult CHGMULT CHGMULT
                        Molecular charge and multiplicity for geometry
  -r [RECIPE ...], --recipe [RECIPE ...]
                        NICS recipe(s) to follow
  -f [{xyz,qchem,g16} ...], --filetype [{xyz,qchem,g16} ...]
                        File types to write augmented molecules to
  -m MODELCHEM, --modelchem MODELCHEM
                        Combination of method/basis to use when preparing input file(s)
  -l LOCATION, --location LOCATION
                        Path to directory within which to write files
  -n NAME, --name NAME  Manually specify output filename
  -e [EXTRACT ...], --extract [EXTRACT ...]
                        Output file(s) from which to extract NICS values
```

#### Example 01: Place NICS probes onto a molecule stored in an XYZ file

For a cyclic molecule with Cartesian coordinates stored in `molecule.xyz`, use the following
command to place NICS(0,+1,-1) probes and write them to a new XYZ file:

```bash
(nicshelper) user@host:~/path/to/repo$ python nics-helper.py -g samples/pbenzyne.xyz -r isosimple -f xyz -l samples/
```


