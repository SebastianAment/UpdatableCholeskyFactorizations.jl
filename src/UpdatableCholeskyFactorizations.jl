module UpdatableCholeskyFactorizations

using LinearAlgebra
using LinearAlgebra: checksquare

export UpdatableCholesky, UCholesky,
       updatable_cholesky, ucholesky,
       updatable_cholesky!, ucholesky!,
       add_column!, addcol!,
       remove_column!, remcol!

# IDEA: have separate package PermutationMatrices.jl,
# where we take care of all the permutation logic
# IDEA: add Inverse to have efficient ternary dot product and multiply
mutable struct UpdatableCholesky{T, UT <: AbstractMatrix{T}} <: Factorization{T}
    U_full::UT # pre-allocated U matrix
    n::Int # size of factorization
    info::Int
end
const UCholesky = UpdatableCholesky

"""
```
updatable_cholesky(A::AbstractMatrix, m::Int = 2size(A, 1); check::Bool = true)
```
Computes an updatable Cholesky factorization starting with the `n by n` matrix `A`,
and pre-allocates an `m x m` matrix for future additions to the matrix (defaults to twice the size `2n`).
"""
function updatable_cholesky(A::AbstractMatrix, m::Int = 2size(A, 1); check::Bool = true)
    n = checksquare(A)
    m ≥ n || throw(DimensionMismatch("$m < size(A, 1) == $n"))
    U_full = zeros(eltype(A), m, m)
    @. U_full[1:n, 1:n] = A
    updatable_cholesky!(U_full, n, check = check)
end
const ucholesky = updatable_cholesky

"""
```
updatable_cholesky!(U_full::AbstractMatrix, n::Int; check::Bool = true)
```
Computes an updatable Cholesky factorization. Assumes that the currently available
data `A` is in the upper left `n x n` submatrix.
"""
function updatable_cholesky!(U_full::AbstractMatrix, n::Int; check::Bool = true)
    m = checksquare(U_full)
    A = @view U_full[1:n, 1:n]
    C = cholesky!(A, check = check)
    @. U_full[1:n, 1:n] = C.U
    UpdatableCholesky(U_full, n, C.info)
end
const ucholesky! = updatable_cholesky!

# constructing
function LinearAlgebra.Cholesky(F::UpdatableCholesky)
    Cholesky(F.U, :U, F.info)
end

function Base.getproperty(F::UpdatableCholesky, s::Symbol)
    if s == :U
        U = @view F.U_full[1:F.n, 1:F.n]
        UpperTriangular(U)
    elseif s == :L
        F.U'
    else
        getfield(F, s)
    end
end
Base.eltype(F::UpdatableCholesky{T}) where T = T
Base.size(F::UpdatableCholesky) = (F.n, F.n)
Base.size(F::UpdatableCholesky, i::Int) = i > 2 ? 1 : size(F)[i]
Base.AbstractMatrix(F::UpdatableCholesky) = Matrix(F)
Base.adjoint(C::UpdatableCholesky) = C # since it is Hermitian
LinearAlgebra.issuccess(C::UpdatableCholesky) = C.info == 0
function Base.Matrix(F::UpdatableCholesky)
    U = F.U
    U'U
end

function Base.copy(F::UpdatableCholesky)
    UpdatableCholesky(copy(F.U_full), F.n, F.info)
end

# linear solves with updatable cholesky
LinearAlgebra.:\(C::UpdatableCholesky, x::AbstractVector) = ldiv!(C, copy(x))
LinearAlgebra.:\(C::UpdatableCholesky, X::AbstractMatrix) = ldiv!(C, copy(X))
LinearAlgebra.:/(X::AbstractMatrix, C::UpdatableCholesky) = rdiv!(copy(X), C)
LinearAlgebra.rdiv!(X::AbstractMatrix, F::UpdatableCholesky) = ldiv!(F, X')'
function LinearAlgebra.ldiv!(F::UpdatableCholesky, x::AbstractVecOrMat)
    U = F.U
    ldiv!(U', x)
    ldiv!(U, x)
    return x
end

"""
```
add_column!(C::UpdatableCholesky, A::AbstractMatrix)
```
appending columns in A - and their corresponding rows A' - to `UpdatableCholesky`
 factorization `C`. Complexity is `O(m³ + n²m)` where `n = size(C, 1)` and `m = size(A, 2)`.
"""
function add_column!(C::UpdatableCholesky, A::AbstractMatrix)
    ensure_enough_memory!(C, A)
    n, m = size(C, 1), size(A, 2)
    new_n = n + m
    if size(A, 1) ≠ new_n
        throw(DimensionMismatch("New columns `A` of $(size(A)) do not have the required dimensionality $new_n to add to UpdatableCholesky of size $n."))
    end
    A12 = @view A[1:n, :]
    A22 = @view A[n+1:end, :]

    S12 = @view C.U_full[1:n, n+1:new_n]
    ldiv!(S12, C.U', A12) # this is squared in the number of old indices and linear in m: O(n²m)

    # S22 = cholesky!(A22 - S12' * S12).U
    S22 = @view C.U_full[n+1:new_n, n+1:new_n]
    @. S22 = A22
    mul!(S22, S12', S12, -1, 1)
    C22 = cholesky!(S22) # this is cubic in the number of new indices only `O(m³)`
    C.n = new_n
    return C
end
function add_column!(C::UpdatableCholesky, a::AbstractVector)
    add_column!(C, reshape(a, :, 1))
end

"""
```
remove_column!(C::UpdatableCholesky, i::Int)
```
updating `UpdatableCholesky` factorization `C` corresponding to removal of the
`i`th row and column of the original matrix. Complexity is `O(n²)`.
"""
function remove_column!(C::UpdatableCholesky, i::Int)
    n = C.n
    if i == n # removing the last row / column is trivial
        C.n -= 1
        return C
    end

    # Create vector for rank update
    t = zeros(eltype(C), n)
    # t = @view C.U[:, i] # IDEA: can we use existing memory? Problem: Base's lowrankupdate! does not accept views.
    @. t[1:i] = 0
    @. t[i+1:n] = conj.(C.U[i, i+1:n])
    # Do rank update
    D = Cholesky(C)
    lowrankupdate!(D, t)

    @. C.U[:, i] = 0 # shift memory to the left, IDEA: could instead update view, but this will be more complex with subsequent additions
    for j in i+1:n
        @. C.U[1:i-1, j-1] = C.U[1:i-1, j]
        @. C.U[i:j-1, j-1] = C.U[i+1:j, j]
    end
    C.n -= 1
    return C
end

const addcol! = add_column!
const remcol! = remove_column!

"""
```
ensure_enough_memory!(C::UpdatableCholesky, A::AbstractMatrix, memory_expansion_factor::Int = 2)
```
makes sure we have enough memory pre-allocated in U_full to accomodate
additions of columns in `A`. `memory_expansion_factor` determines how much memory
is requested to accomodate future column additions, if we don't have enough space
to add `A`.
"""
function ensure_enough_memory!(C::UpdatableCholesky, A::AbstractMatrix, memory_expansion_factor::Int = 2)
    if size(C.U_full, 1) < C.n + size(A, 2) # create a bigger U_full, and store it in C.U_full
        n = C.n
        m = memory_expansion_factor * n # new size
        U_full = zeros(eltype(C), m, m)
        @. U_full[1:n, 1:n] = C.U
        C.U_full = U_full
    end
    return C
end

end # module
