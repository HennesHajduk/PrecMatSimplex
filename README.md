# Sparsified discrete gradient operators

This repository contains codes to compute Bernstein finite element operators on the reference simplex.
These matrices represent discrete gradient operators that are sparsified such that only entries
corresponding to closest neighbors and diagonal entries are nonzero.

## Old sparsification

This sparsification is achieved by multiplying consistent dense discrete gradient operators from the
left by the reference element lumped mass matrix times the inverse of the reference element consistent
mass matrix. Due to ill-conditioning of high-order Bernstein matrices, these preconditioned matrices
should not be computed by inversion, but rather using information of the sparsified element stencil
from the proof of the fact that preconditioning results in sparsification, see [Kuzmin and Quezada de
Luna (2020) JCP, Appendix](https://www.sciencedirect.com/science/article/pii/S0021999120301856).

The codes in the directory ***non-skew-symmetric-triangle*** compute the version proposed in that
paper for triangles. It represents the state of this GitHub repository upon publication of the article
[Hajduk (2021) CAMWA](https://www.sciencedirect.com/science/article/abs/pii/S0898122121000547).

The routine *computePrecMatSimplex* in the main directory generalizes this code to simplices in one,
two, and three dimensions. The routine *computePrecMat1D* uses the approach first proposed in [Lohmann
et al. (2017) JCP, Eq. (B.4)](https://www.sciencedirect.com/science/article/pii/S0021999117303388) for
the 1D case. Note that the setting as required therein is slightly different from ours (their operator
is equal to our matrix if v=-1 in that reference). Using tensor products, the multidimensional version
of this matrix is computed by the routine *computePrecMatBox*.

### Citation

If you use this code in your work, please cite
[Hajduk (2021) CAMWA](https://www.sciencedirect.com/science/article/abs/pii/S0898122121000547).

You may also want to cite [Lohmann et al. (2017)
JCP](https://www.sciencedirect.com/science/article/pii/S0021999117303388)
if you use 1D elements, quadrilaterals, or hexahedra, and/or [Kuzmin and Quezada de Luna (2020)
JCP](https://www.sciencedirect.com/science/article/pii/S0021999120301856) for simplices or in general.

## New sparsification

The old sparsification approach does not produce operators with skew-symmetric off-diagonal entries.
This property may, however, be desirable for multiple reasons. Thus, a new version providing skew
symmetry has been proposed recently.
The routines *computeSkewDiscreteGradientGeometry* compute these operators for all supported element
types, where *Geometry* is either *1D*, *Box*, *Prism*, or *Simplex*.

### Citation

If you use this code in your work, please cite **Hajduk** (2024) ***Improvements of algebraic flux-
correction schemes based on Bernstein finite elements***, in review: Journal of Numerical Mathematics.
