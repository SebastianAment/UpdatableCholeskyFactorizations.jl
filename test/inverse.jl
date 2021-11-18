module TestInverseCholeskyDot
using LinearAlgebra
using UpdatableCholeskyFactorizations
using LazyInverses
using Test

# testing particularly efficient algebraic operations with inverse(Cholesky)
@testset "inverse" begin
    k = 16
    element_types = (Float32, Float64, ComplexF32, ComplexF64)
    for elty in element_types
        for n in (16, 1024) # to trigger sequential and parallel code
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
            @test diag(invC) ≈ diag(inv(A))
        end
    end
end

end # TestInverseCholeskyDot
