# adding Inverse to have efficient ternary dot product and multiply
# allows for especially efficient Ma
import LinearAlgebra: dot, *, diag
function dot(x::AbstractVector, A::Inverse{<:Any, <:UpdatableCholesky}, y::AbstractVector)
	xp = copy(x)
	yp = x ≡ y ? xp : copy(y)
	dot!!(xp, A, yp)
end
# IDEA: could have parallel_threshold argument here
import LazyInverses: dot!!, inverse_cholesky_dot!!
function dot!!(x::AbstractVector, A::Inverse{<:Any, <:UpdatableCholesky}, y::AbstractVector)
	C = A.parent
	inverse_cholesky_dot!!(x, C, y)
end

using LazyInverses: inverse_cholesky_diag, inverse_cholesky_mul, inverse_cholesky_diag_mul
function diag(A::Inverse{<:Any, <:UpdatableCholesky})
	inverse_cholesky_diag(A.parent)
end

import LazyInverses: diag_mul
function diag_mul(X, A::Inverse{<:Any, <:UpdatableCholesky}, Y)
	C = A.parent
	inverse_cholesky_diag_mul(X, C, Y)
end

function *(X, A::Inverse{<:Any, <:UpdatableCholesky}, Y)
	inverse_cholesky_mul(X, A.parent, Y)
end
