# adding Inverse to have efficient ternary dot product and multiply
# allows for especially efficient Ma
import LinearAlgebra: dot
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
