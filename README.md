# UpdatableCholeskyFactorizations.jl
[![CI](https://github.com/SebastianAment/UpdatableCholeskyFactorizations.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/SebastianAment/UpdatableCholeskyFactorizations.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/SebastianAment/UpdatableCholeskyFactorizations.jl/branch/main/graph/badge.svg?token=8L49I28I73)](https://codecov.io/gh/SebastianAment/UpdatableCholeskyFactorizations.jl)

This package contains implementations of efficient representations and updating algorithms for Cholesky factorizations.
Compared to constructing Cholesky factorizations from scratch, this package can lead to significant speed improvements, see the benchmarks below.

## Installation
To install the package, type `]` and subsequently, `add UpdatableCholeskyFactorizations` in the Julia REPL.

## Basic Usage
The package exports the `UpdatableCholesky` type, which can also be referenced as `UCholesky`. 
A structure of this type can be computed using the function
```julia
updatable_cholesky(A::AbstractMatrix, m::Int = 2size(A, 1); check::Bool = true),
```
which computes an updatable Cholesky factorization starting with the `n x n` matrix `A`,
and pre-allocates an `m x m` matrix for future additions to the matrix (defaults to twice the size `2n`).
Given an `UpdatableCholesky` structure,
we can use
```julia
add_column!(C::UpdatableCholesky, A::AbstractMatrix)
```
to append columns in A - and their corresponding rows A' - to `UpdatableCholesky` factorization `C`. Complexity is `O(m³ + n²m)` where `n = size(C, 1)` and `m = size(A, 2)`.
To remove a column from `C`,
```julia
remove_column!(C::UpdatableCholesky, i::Int)
```
updates the `UpdatableCholesky` factorization `C` corresponding to removal of the
`i`th row and column of the original matrix. Complexity is `O(n²)`, where `n = size(C, 1)`.
See the following example.

## Example
```julia
using LinearAlgebra
using UpdatableCholeskyFactorizations

# setting up 
n = 128
m = 16
A_full = randn(2n, 2n)
A_full = A_full'A_full
A = @view A_full[1:n, 1:n]
a = @view A_full[1:n+1, n+1]
Bm = @view A_full[1:n+m, n+1:n+m] # m columns

C = updatable_cholesky(A) # computing a factorization of A

add_column!(C, a) # appends the column a to the factorization, C is now a factorization of A_full[1:n+1, 1:n+1]

remove_column!(C, n+1) # removes the (n+1)st column from the factorization, C is now a factorization of A_full[1:n, 1:n]

add_column!(C, Bm) # appends the m columns in Bm to the factorization, C is now a factorization of A_full[1:n+m, 1:n+m]

```

## Efficiency

The following benchmarks highlight the performance benefits that updating a Cholesky factorization can have over constructing a new factorization from scratch.
First, we report the performance of `LinearAlgebra`'s `cholesky` on two matrix sizes:
```julia
128 by 128 cholesky
BenchmarkTools.TrialEstimate: 
  time:             85.072 μs
  gctime:           0.000 ns (0.00%)
  memory:           128.12 KiB
  allocs:           4
  
256 by 256 cholesky
BenchmarkTools.TrialEstimate: 
  time:             315.872 μs
  gctime:           0.000 ns (0.00%)
  memory:           512.12 KiB
  allocs:           4
```
Now, compare this to this package's `updatable_cholesky` and `add_column!` functions:
```julia
128 by 128 updatable_cholesky
BenchmarkTools.TrialEstimate: 
  time:             169.834 μs
  gctime:           0.000 ns (0.00%)
  memory:           512.44 KiB
  allocs:           8

adding vector to 128 by 128 updatable_cholesky
BenchmarkTools.TrialEstimate: 
  time:             4.419 μs
  gctime:           0.000 ns (0.00%)
  memory:           208 bytes
  allocs:           3

adding 16 vectors to 128 by 128 updatable_cholesky
BenchmarkTools.TrialEstimate: 
  time:             23.295 μs
  gctime:           0.000 ns (0.00%)
  memory:           208 bytes
  allocs:           3
```
The construction of the `UpdatableCholesky` structure tends to be approximately two times slower than `cholesky` for a given size,
and by default, consumes as much memory as a cholesky factorization of twice the size, in order to accomodate future column additions without requiring additional memory allocations.
However, subsequently updating a factorization is then much faster than computing another factorization from scratch and does not require additional allocations.
In the benchmarks above, adding a vector to the 128 by 128 factorization is ~**20 times faster**,
and adding 16 vectors is ~**4 times faster** than calling `cholesky`.
The benchmarks were computed on a 2017 MacBook Pro with an i7 dual core and 16 GB of RAM.

## Future Work
Currenlty, the package only supports appending columns. Future work could add support for adding columns at arbitrary indices.


