## Class definitions

## Accepted input representations for the matrix-signature methods.
##
## `Matrix` is the virtual base class of the Matrix package, so this one
## member covers every representation it defines -- dgCMatrix, dgRMatrix,
## dgTMatrix, the dense dgeMatrix, the logical lgCMatrix -- instead of
## just the column-compressed one.
##
## Not covered: `SVT_SparseMatrix` from the SparseArray package, which is
## not a `Matrix` subclass. Supporting it would add a dependency, so it
## is deliberately out of scope here.
setClassUnion(name = "AnyMatrix", members = c("matrix", "Matrix"))
