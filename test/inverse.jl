module TestInverseCholeskyDot
using LinearAlgebra
using UpdatableCholeskyFactorizations
using LazyInverses
using UpdatableCholeskyFactorizations: diag_mul
using Test

# testing particularly efficient algebraic operations with inverse(Cholesky)
@testset "inverse" begin
    k = 16
    element_types = (Float32, Float64, ComplexF32, ComplexF64)
    for elty in element_types
        @testset "eltype $elty" begin
            for n in (16) #, 1024) # to trigger sequential and TODO: parallel code
                @testset "size $n" begin
                    x, y = rand(elty, n), rand(elty, n)
                    X, Y = rand(elty, k, n), rand(elty, n, k)
                    A = randn(elty, n, n)
                    A = A'A + I # for numerical stability, especially for low-precision floats
                    # with UpdatableCholesky
                    C = updatable_cholesky(A)
                    invC = Inverse(C)
                    @test dot(x, invC, x) ≈ dot(x, C\x)
                    @test dot(x, invC, y) ≈ dot(x, C\y)
                    @test *(X, invC, X') ≈ *(X, C\X')
                    @test *(X, invC, Y) ≈ *(X, C\Y)
                    @test diag_mul(X, invC, X') ≈ diag(*(X, invC, X'))
                    @test diag_mul(X, invC, Y) ≈ diag(*(X, invC, Y))
                    @test diag(invC) ≈ diag(inv(A))
                end # testset n
            end # loop over n
        end # testset eltype
    end # loop over eltype
end

end # TestInverseCholeskyDot
