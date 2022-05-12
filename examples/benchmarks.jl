using UpdatableCholeskyFactorizations
using LinearAlgebra
using BenchmarkTools

elty = Float64
n = 128
m = 16
A_full = randn(elty, 2n, 2n)
A_full = A_full'A_full
A = A_full[1:n, 1:n]
a = A_full[1:n+1, n+1]
B = A_full[1:2n, n+1:2n]
Bm = A_full[1:n+m, n+1:n+m] # adding m columns
C = updatable_cholesky(A)

set = @benchmarkset "updatable cholesky" begin # for n in 1:16 IDEA: could have loop over different sizes
    @case "$n by $n cholesky" cholesky($A)
    @case "$(2n) by $(2n) cholesky" cholesky($A_full)
    @case "$n by $n updatable_cholesky" updatable_cholesky($A)
    @case "adding vector to $n by $n UpdatableCholesky" add_column!(D, $a) setup=(D=copy($C))
    @case "adding $m vectors to $n by $n UpdatableCholesky" add_column!(D, $Bm) setup=(D=copy($C))
    @case "adding $n x $n matrix to $n by $n UpdatableCholesky" add_column!(D, $B) setup=(D=copy($C))
end

result = run(set)

r = minimum(result["updatable cholesky"])
# times = Vector{Float64}(undef, length(r))
# titles = Vector{String}(undef, length(r))
# i = 1
for (title, trial) in r.data
    println()
    println(title)
    display(trial)
    # times[i] = trial.time
    # titles[i] = title
    # i += 1
end
