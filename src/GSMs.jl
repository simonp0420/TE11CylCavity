module GSMs
export GSM, cascade, cascade!, propagate!

using LinearAlgebra
using StaticArrays: SMatrix

"""
    GSM{T<:Number, M<:AbstractMatrix{T}} <: Any
A struct for storing the partitions of a GSM (Generalized Scattering Matrix)

## Fields:
- `S11`, `S12`, `S21`, `S22`: The partitions of the GSM.  Each is of type `AbstractMatrix{T}`
"""
struct GSM{T<:Number, M<:AbstractMatrix{T}}
    s11::M
    s12::M
    s21::M
    s22::M
end

"""
    GSM(n1::Int, n2::Int)

Convenience constructor.
## Arguments
- `n1`, `n2`: These specify the dimensions of the GSM partition matrices.
  `s11` will have size `(n1,n1)`, `s21` will have size `(n2,n1)`, etc. The partition
  matrices will be of type `Matrix{ComplexF64}`.  `s11` and `s12` will be initialized
  to all zeros, while `s12` and `s21` will be initialized to identity matrices.
"""
function GSM(n1::Int, n2::Int)
    gsm = GSM(zeros(ComplexF64, n1, n1), zeros(ComplexF64, n1, n2),
        zeros(ComplexF64, n2, n1), zeros(ComplexF64, n2, n2))
    gsm.s12[diagind(gsm.s12)] .= one(eltype(gsm.s12))
    gsm.s21[diagind(gsm.s21)] .= one(eltype(gsm.s21))
    gsm
end

"""
    GSM(s11::Number, s12::Number, S21::Number, s22::Number)

Convenience constructor for two-port network.
## Arguments
- `s11`, `s12`, `s21`, `s22`: The scalar entries in the scattering matrix. They can be entered either
  as positional or keyword arguments.
"""
function GSM(s11::Number, s12::Number, s21::Number, s22::Number)
    s11mat = SMatrix{1, 1, ComplexF64, 1}(s11)
    s12mat = SMatrix{1, 1, ComplexF64, 1}(s12)
    s21mat = SMatrix{1, 1, ComplexF64, 1}(s21)
    s22mat = SMatrix{1, 1, ComplexF64, 1}(s22)
    gsm = GSM(s11mat, s12mat, s21mat, s22mat)
    return gsm
end

GSM(; s11::Number, s12::Number, s21::Number, s22::Number) = GSM(s11, s12, s21, s22)


@inline function Base.getindex(gsm::GSM, i, j)
    (i, j) == (1, 1) && (return gsm.s11)
    (i, j) == (1, 2) && (return gsm.s12)
    (i, j) == (2, 1) && (return gsm.s21)
    (i, j) == (2, 2) && (return gsm.s22)
    throw(BoundsError(gsm, (i, j)))
end

Base.size(::GSM) = (2, 2)
function Base.size(::GSM, k)
    (k == 1 || k == 2) && (return 2)
    k > 2 && (return 1)
    error("arraysize: dimension out of range")
end


"""
    cascade(a::GSM, b::GSM) -> c::GSM

Cascade a pair of generalized scattering matrices (GSMs).

## Input arguments

- `a`, `b`:  variables containing the input GSMs that are to be cascaded.  The matrices 
must be conformable, i.e., `n2a == n1b`, where `n2a` is the number of modes in Region
2 for GSM `a`, and `n1b` is the number of modes in Region 1 of GSM `b`.
"""
function cascade(a::GSM, b::GSM)
    n2a = size(a[2,2], 2)
    n1b = size(b[1,1], 1)
    n1b ≠ n2a && error("Non-conformable arrays")
    # Equation (3.35) of the theory documentation:
    gprod1 = (I - a[2,2] * b[1,1]) \ a[2,1]
    s21 = b[2,1] * gprod1
    s11 = a[1,1] + (a[1,2] * b[1,1] * gprod1)
    gprod2 = (I - b[1,1] * a[2,2]) \ b[1,2]
    s12 = a[1,2] * gprod2
    s22 = b[2,2] + (b[2,1] * a[2,2] * gprod2)
    return GSM(s11, s12, s21, s22)
end



end # module
