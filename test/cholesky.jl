module TestUpdatableCholeskyFactorizations
using UpdatableCholeskyFactorizations
using LinearAlgebra
using Test

@testset "updating Cholesky" begin
    element_types = (Float32, Float64, ComplexF32, ComplexF64)
    for elty in element_types
        n = 8
        A_full = randn(elty, 2n, 2n)
        A_full = A_full'A_full
        A = @view A_full[1:n, 1:n]

        UC = updatable_cholesky(A)
        @test UC isa UpdatableCholesky{elty}
        C = cholesky(A) # compare against regular cholesky
        @test UC.U ≈ C.U
        @test UC.L ≈ C.L
        @test UC.info == C.info
        @test size(UC) == size(A)
        for i in 1:3
            @test size(UC, i) == size(A, i)
        end

        # conversion to Cholesky type
        C2 = Cholesky(UC)
        @test C2 isa Cholesky{elty}
        @test C2.U ≈ C.U
        @test C2.L ≈ C.L
        @test C2.info == C.info
        @test size(C2) == size(C)

        # adding rows and columns
        for i in 1:n
            a = @view A_full[1:n+i, n+i]
            add_column!(UC, a)
            @test Matrix(UC) ≈ @view A_full[1:n+i, 1:n+i]
            @test size(UC) == (n+i, n+i)
        end

        # removing rows and columns until nothing is left
        for i in reverse(1:2n)
            remove_column!(UC, i)
            @test Matrix(UC) ≈ @view A_full[1:i-1, 1:i-1]
            @test size(UC) == (i-1, i-1)
        end

        # add columns in again
        for i in 1:n
            a = @view A_full[1:i, i]
            add_column!(UC, a)
            @test Matrix(UC) ≈ @view A_full[1:i, 1:i]
            @test size(UC) == (i, i)
        end

        # and remove columns, this time in different order
        for i in 1:n
            remove_column!(UC, 1) # keep removing from the front
            @test Matrix(UC) ≈ @view A_full[1+i:n, 1+i:n]
            @test size(UC) == (n-i, n-i)
        end
    end
end

end
