module UpdatableCholeskyFactorizations

using LinearAlgebra
using LinearAlgebra: checksquare

export UpdatableCholesky, UCholesky,
       updatable_cholesky, ucholesky,
       updatable_cholesky!, ucholesky!,
       add_column!, addcol!,
       remove_column!, remcol!

using LazyInverses

# IDEA: have separate package PermutationMatrices.jl,
# where we take care of all the permutation logic

include("cholesky.jl")
include("inverse.jl")

end # module
