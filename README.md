_Supplementary Material for_ Variations on the Bergman Cyclization Theme: Electrocyclizations of Penta-, Hepta-, and Octa-Diynes
===================================================================================================

Welcome to the supplementary repository for the 5-7-8 project! Here, we make
available the software that was written to facilitate the automatic generation,
performance, and analysis of nucleus-independent chemical shift computations
during the course of this project. 

Enjoy!

~~ The Authors

**Please use the following citation when referencing this article and/or using
this software:**


## Motivation

As discussed in both the manuscript and Supplemental Information for this work,
while NICS is a very useful index to assess the aromaticity of a molecular
system, it is unfortunately sensitive to the location at which the probes are
placed. Furthermore, for the non-symmetric, non-planar cyclic molecules
examined here, conventional NICS probe locations are not well-defined.
Therefore, we have developed a robust, general approach to place NICS probes
for even non-symmetric and non-planar molecules utilizing the principal moment
of inertia of the non-mass--weighted ring atoms. In addition to being
applicable to a more diverse set of molecules, this approach enjoys
self-consistency with the previous literature leveraging NICS because it
reproduces the conventional probe locations for planar, symmetric molecules by
construction.

## Theory

Do we want an overview of the theory here? Or is it just sufficient to refer to
Section I.G of the SI where I derive it?

## Using the `nics-helper.py` Script

### Installation

To use the helper script `nics-helper.py`, first clone this repository to
your local machine before proceeding to install the required dependencies.

It is recommended that you perform this installation in a new Conda environment
to protect against version clashes with other Python packages you may have
installed on your machine. To do this, first download and install
[Anaconda](https://www.anaconda.com/products/distribution) or
[Miniconda](https://docs.conda.io/en/latest/miniconda.html), before executing
the following from the command line:

```bash
user@host:~$ conda create -n nicshelper python=3.9 qcelemental networkx -c conda-forge
```

This command will create a new Conda environment named `nicshelper`, into which
all the dependencies for the `nics-helper.py` script will be installed. To then
use the script, first activate the Conda environment you just created:

```bash
user@host:~$ conda activate nicshelper
```

And voila! Now you're ready to use the script.

### Examples

The `nics-helper.py` script can be used either by importing it as a Python
module into another script or by calling it directly from the command line.
Below you'll find some typical use cases and how to call the script from the
command line to perform those actions (using representative sample files 
available in the `samples/` directory).

#### Example 00: Getting a help summary from the script

Once you've activated the `nicshelper` Conda environment, run the following
from the command line to get a summary of the command-line arguments to
`nics-helper.py`:

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

For the cyclic diradical molecule _p_-benzyne with Cartesian coordinates stored
in `samples/pbenzyne.xyz`, use the following command to place NICS(0,+1,-1)
probes and write them to a new XYZ file for visualization:

```bash
(nicshelper) user@host:~/path/to/repo$ python nics-helper.py -g samples/pbenzyne.xyz -r isosimple -f xyz -l samples/
```

Executing this command in the top-level repository directory of your clone will
generate the file `samples/pbenzyne-nicsviz.xyz`, with the following contents:

```
13
0 1 C6H4He3
C                     1.207143108871    -0.708245481697     0.000000000000
H                     2.146410188946    -1.229570401074     0.000000000000
C                     0.000000000000    -1.340286751545     0.000000000000
C                    -1.207143108871    -0.708245481697     0.000000000000
H                    -2.146410188946    -1.229570401074     0.000000000000
C                     1.207143108871     0.708245481697     0.000000000000
H                     2.146410188946     1.229570401074     0.000000000000
C                     0.000000000000     1.340286751545     0.000000000000
C                    -1.207143108871     0.708245481697     0.000000000000
H                    -2.146410183654     1.229570401074     0.000000000000
He                    0.000000000000     0.000000000000    -1.000000002404
He                    0.000000000000     0.000000000000     0.000000000000
He                    0.000000000000     0.000000000000     1.000000002404
```

Here, the NICS probes are represented as He atoms to facilitate their
visualization in a GUI (e.g., [IQMol](http://iqmol.org/downloads.html)).

>*Note*: Probes are not internally labeled (i.e., as corresponding to NICS(-1),
>NICS(0), NICS(+1), etc.probe); however they are _always_ added/printed in
>ascending order. So, in the XYZ file above, NICS(-1) is printed first, above
>NICS(0), and finally above NICS(+1).

#### Example 02: Generate a Gaussian16 input file for NICS computation on _p_-benzyne

If you want to do more than just _place_ NICS probes, never fear!
`nics-helper.py` can generate both Gaussian16 and Q-Chem input files for NICS
computations. To demonstrate this capability, let's generate the Gaussian input
file for our cyclic diradical, _p_-benzyne:

```bash
(nicshelper) user@host:~/path/to/repo$ python nics-helper.py -g samples/pbenzyne.xyz -r isosimple -f g16 -c 0 1 -m BS-UB3LYP/6-311++G** -l samples/
```

This will generate the following input file, written to `samples/pbenzyne_geom-BSUB3LYP_6-311ppGss_nics.com`:

```
%chk=pbenzyne_geom-BSUB3LYP_6-311ppGss_nics.chk
%nproc=2
#UB3LYP/6-311++G** guess=(mix,always) scf=qc nmr=giao

 C6H4X3: Broken-symmetry requested; mixing HOMO/LUMO in SCF guess before NICS job

0 1
C                     2.281169870000    -1.338389990000     0.000000000000
H                     4.056127410000    -2.323551310000     0.000000000000
C                     0.000000000000    -2.532774890000     0.000000000000
C                    -2.281169870000    -1.338389990000     0.000000000000
H                    -4.056127410000    -2.323551310000     0.000000000000
C                     2.281169870000     1.338389990000     0.000000000000
H                     4.056127410000     2.323551310000     0.000000000000
C                     0.000000000000     2.532774890000     0.000000000000
C                    -2.281169870000     1.338389990000     0.000000000000
H                    -4.056127400000     2.323551310000     0.000000000000
Bq                    0.000000000000     0.000000000000    -1.889726130000
Bq                    0.000000000000     0.000000000000     0.000000000000
Bq                    0.000000000000     0.000000000000     1.889726130000

```

In this example, because our molecule is a diradical, we have broken the spin
and spatial symmetry by mixing the HOMO and LUMO throughout the SCF procedure;
this is accomplished with the `guess=(mix,always)` option.

>***Warning:*** If your molecule is a diradical and you are targeting the open-shell
>singlet electronic configuration for your NICS job:
>    * You _must_ use Gaussian to perform the NICS computations, and
>    * You _must_ request the broken symmetry approach by prepending `BS-` onto
>    the level of theory (combination of method/basis set) passed to the
>    `--modelchem` command-line option of `nics-helper.py.
>
>Otherwise, the computation will produce garbage results!

>Note: The `scf=qc` option changes the SCF convergence algorithm to one which is
>more robust than the default, and therefore more likely to converge (albeit
>slightly more slowly).

#### Example 03: Extracting NICS data from Gaussian output

After performing a NICS computation, `nics-helper.py` can also make extracting NICS
data from the output file a snap!

To extract NICS data from the output of a Gaussian computation
`samples/pbenzyne_geom-BSUB3LYP_6-311ppGss_nics.log`, pass the `-e` option to
`nics-helper.py`:

```bash
(nicshelper) user@host:~/path/to/repo$ python nics-helper.py -e samples/pbenzyne_geom-BSUB3LYP_6-311ppGss_nics.log
```

This will write a new XYZ file alongside the output, named
`samples/pbenzyne_geom-BSUB3LYP_6-311ppGss_nicsdata.xyz`, with the following
contents:

```
13
Data extracted from output file: samples/pbenzyne_geom-BSUB3LYP_6-311ppGss_nics.log
C     2.2811700000 -1.3383900000  0.0000000000
H     4.0561270000 -2.3235510000  0.0000000000
C     0.0000000000 -2.5327750000  0.0000000000
C    -2.2811700000 -1.3383900000  0.0000000000
H    -4.0561270000 -2.3235510000  0.0000000000
C     2.2811700000  1.3383900000  0.0000000000
H     4.0561270000  2.3235510000  0.0000000000
C     0.0000000000  2.5327750000  0.0000000000
C    -2.2811700000  1.3383900000  0.0000000000
H    -4.0561270000  2.3235510000  0.0000000000
X     0.0000000000  0.0000000000 -1.8897260000  0.1592000000
X     0.0000000000  0.0000000000  0.0000000000  2.5557000000
X     0.0000000000  0.0000000000  1.8897260000  0.1592000000
```

Here, the fourth column contains the isotropic deshielding values experienced
at each NICS probe (that is, the negative of the isotropic shielding values
reported in the output file!). Additionally, the probe symbol is `X` because
the fourth column breaks XYZ file format sufficiently to nullify opening it in
a GUI.  These deshielding values indicate the relative aromaticity of the
molecule:
* `deshielding` < -5 : Aromatic
* -5 <= `deshielding` <= 5 : Non-aromatic
* 5 < `deshielding` : Anti-aromatic


