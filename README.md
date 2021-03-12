# PrecMatSimplex
Preconditioned gradient matrix on the reference simplex

This code computes the preconditioned gradient matrix on the reference simplex according to

https://www.sciencedirect.com/science/article/abs/pii/S0898122121000547

section 3, particularly eq. (3.2) and the appendices.

The preconditioning is achieved by multiplying the gradient element matrices from the left with the reference element lumped mass matrix times the inverse of the reference element consistent mass matrix.

If the finite element functions are represented in the Bernstein basis, this preconditioning produces sparsified element matrices with non-zero entries corresponding only to nearest neighbor nodes. Due to the ill-conditioning of high-order Bernstein mass matrices, these preconditioned matrices should not be computed by inversion, but rather using information of the sparsified element stencil from the proof of the fact that precondioning results in sparsification.

See also the respective appendices in 

https://www.sciencedirect.com/science/article/pii/S0021999117303388

https://www.sciencedirect.com/science/article/pii/S0021999120301856
