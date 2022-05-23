`iwpert`: Perturb Structure Along Vibrational Mode(s)
=====================================================

### Overview

One key component of molecular geometry optimization is confirming the
identified structures are indeed minima (for stable species) or first-order
saddle points (for transition states) on the potential energy surface. This is
done via frequency analysis, whereby the vibrational frequencies and normal
modes are computed by diagonalizing the molecular Hessian matrix.  For minima
and transition states, this process must return all real vibrational
frequencies or exactly one imaginary frequency, respectively. But what if there
are more imaginary modes than are desired? Or what if the imaginary mode
identified by a TS search is not the correct one?

Both of these situations may be remedied by projecting out any undesired
imaginary mode(s) by adding a scalar multiple of the imaginary mode vector(s)
to the Cartesian geometry of the molecule. This can be straightforwardly
accomplished using a spreadsheet software, however when many such operations
are necessary, this manual approach becomes rapidly tedious and is prone
to copy-paste errors throughout.

Instead, we have developed a Python script to automatically construct a
perturbed molecular geometry, projected along the desired vibrational mode(s)
by a requested scalar factor, starting simply from the output file of a
frequency computation performed by Q-Chem.

### Using `iwpert.py`



