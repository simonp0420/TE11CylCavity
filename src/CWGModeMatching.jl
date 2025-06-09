"""
    CWGModeMatching

Mode matching module for circular waveguides with m = 1 azimuthal variation.
"""
module CWGModeMatching

export CWG, setup_modes!, cascade


include("J1Jp_roots.jl")

using SpecialFunctions: besselj0, besselj1, besselj
using LinearAlgebra: I, \, transpose, lu!, ldiv!, mul!

include("Constants.jl")
using .Constants: c₀, η₀, μ₀

include("GSMs.jl")
using .GSMs: GSMs, GSM
import .GSMs: cascade

@enum TETM::Bool TE=true TM=false

"""
    Mode
A struct describing a m=1 mode.

## Fields
- `n::Int`: Radial variation mode number.
- `p::TETM`: Mode type (`TE` or `TM`).
- `kcoa::Float64`: Cutoff wavenumber multiplied by waveguide radius.
- `γ::ComplexF64`: Attenuation constant [np/inch]. Lies in first quadrand of complex plane.
- `Z::ComplexF64`: Mode impedance [Ohm].

## Constructor
Constructors for `Mode` can use either keyword or positional arguments. `γ` and `Z` are both optional and will 
default to zero.
"""
@kwdef struct Mode
    n::Int  # radial variation mode index
    p::TETM # Mode type (TE or TM)
    kcoa::Float64 # Cutoff wavenumber * waveguide radius [radian]
    γ::ComplexF64 = zero(ComplexF64) # complex mode attenuation constant [np/inch]
    Z::ComplexF64 = zero(ComplexF64)# complex mode impedance normalized to η₀
end

"""
    CWG
A struct describing a section of metallic circular waveguide.

## Fields
- `a::Float64`: Inner radius [inch]
- `l::Float64`: Section length [inch]
- `ϵᵣ`::Float64`: Dielectric constant of the medium filling the guide.
- `tanδ`::Float64`: Loss tangent of the medium filling the guide.
- `σ::Float64`: Wall conductivity [S/m]
- `modes::Vector{Mode}`: List of modes considered in this section of guide.

## Constructor
Constructors for `Mode` can use either keyword or positional arguments. `σ` and `modes` are both optional and will 
default to `Int` and `Mode[]`, respectively.
"""
@kwdef struct CWG
    a::Float64  # Radius [inch]
    l::Float64  # Length [inch]
    ϵᵣ::Float64 = 1.0    # Relative permittivity of dielectric
    tanδ::Float64  = 0.0 # Loss tangent of dielectric
    σ::Float64 = Inf  # Wall conductivity [S/m]
    modes::Vector{Mode} = Mode[]
end


"""
    mysqrt(x)

Same as `sqrt` unless `sqrt(x)` is pure negative imaginary in which
case it returns `-sqrt(x)` (i.e., positive pure imaginary).
"""
mysqrt(x) = sqrt(x)
function mysqrt(z::Complex)
    ans = sqrt(z)
    return iszero(real(ans)) && imag(ans) < 0 ? -ans : ans
end


J1p(x) = (besselj0(x) - besselj(2, x)) / 2 # J1'(x)
J1xox(x) = besselj0(x) - J1p(x) # J1(x)/x


const chi = first(J1J1p_roots())  # Obtain the Bessel function J1 and J1' roots
const chip = last(J1J1p_roots())  # Obtain the Bessel function J1 and J1' roots

# Evaluate bessel functions J1(chip) and J1'(chi), along with J0 and J2
const J1chip = besselj1.(chip)
const J0chi = besselj0.(chi)
const J0chip = besselj0.(chip)
const J2chi = besselj.(2, chi)

"""
    Iminus!(M)  Replaces matrix `M` with `I-M``.
"""
function Iminus!(M::Matrix)
    m = minimum(size(M))
    @inbounds @simd ivdep for i in eachindex(M)
        M[i] = -M[i]
    end
    @inbounds @simd ivdep for i in 1:m
        M[i, i] += one(eltype(M))
    end
    return M
end

"""
    twoIminus!(M)  replaces matrix `M` with `2*I-M`
"""
function twoIminus!(M::Matrix)
    m = minimum(size(M))
    @inbounds @simd ivdep for i in eachindex(M)
        M[i] = -M[i]
    end
    two = 2 * one(eltype(M))
    @inbounds @simd ivdep for i in 1:m
        M[i, i] += 2
    end
    return M
end

"""
    minusI!(M)  replaces matrix `M` with `M-I`
"""
function minusI!(M::Matrix)
    m = minimum(size(M))
    @inbounds @simd ivdep for i in 1:m
        M[i, i] -= one(eltype(M))
    end
    return M
end

"""
    twoIplus!(M) replaces matrix `M` with `2*I + M`
"""
function twoIplus!(M::Matrix)
    m = minimum(size(M))
    two = 2 * one(eltype(M))
    @inbounds @simd ivdep for i in 1:m
        M[i, i] += 2
    end
    return M
end

"""
    σ2Rs(σmks, fghz) -> Rs

Calculate surface resistance from conductivity and frequency.
## Input Arguments
- `σmks`: Metal bulk conductivity in MKS units [Siemen/meter]
- `fghz`: Frequency in GHz
## Return Value
- `Rs`: Metal surface resistance in Ohms
"""
function σ2Rs(σmks, fghz)
    isinf(σmks) && return 0.0
    ω = 2π * fghz * 1e9 # [radian/s]
    Rs = sqrt(ω * μ₀ / 2σmks) # [Ohm]
    return Rs
end



"""
    setup_modes!(c::CWG, fghz::Float64, nmodes::Int=length(c.modes))

Set up the modes for a uniform circular waveguide.

## Arguments
- `c`: A `CWG` instance.  If `c.modes` is empty, then it will be appended to `length(n)`.  If it is
  already allocated, then its values will be replaced with updated values of `γ` and `Z` corresponding to the 
  new frequency `fghz`.
- `fghz`: The frequency [GHz]
- `nmodes`: (optional) The (even) number of modes to append to `c.modes` if it is empty.  If `c.modes` is nonempty, then
  `nmodes` must be equal to `length(c.modes)`.
"""
function setup_modes!(c::CWG, fghz::Float64, nmodes::Int=length(c.modes))

    if isempty(c.modes) 
        nmodes ≤ 0 && throw(ArgumentError("nmodes must be > 0"))
    else
        nmodes == length(c.modes) || throw(ArgumentError("nmodes not equal to number of existing modes in c"))
    end
    iseven(nmodes) || throw(ArgumentError("nmodes must be even"))

    nmodeso2 = nmodes ÷ 2 # Number of TE modes (= number of TM modes)
    λ₀ = c₀ / fghz
    k₀ = 2π / λ₀
    rootϵ = mysqrt(c.ϵᵣ * complex(1.0, -c.tanδ)) 
    k = k₀ * rootϵ
    k² = k^2
    η = η₀ / rootϵ 
    Rs = σ2Rs(c.σ, fghz)

    if isempty(c.modes)
        # Initialize modes
        for n in 1:nmodeso2,  (p, kcoa) in ((TE, chip[n]), (TM, chi[n]))
            mode = Mode(; n, p, kcoa)
            push!(c.modes, mode)
        end
    end

    # Update γ and Z
    q = 0 # Initialize mode counter
    for n in 1:nmodeso2, p in (TE, TM)
        q += 1
        mode = c.modes[q]
        mode.p == p || error("p mismatch for mode $q!")
        kco = mode.kcoa / c.a
        γ = mysqrt(kco^2 - k²)
        ratio = (kco / real(k))^2
        if !iszero(Rs)
            α = Rs / (c.a * real(η)) / mysqrt(1 - ratio)  # attenuation due to metal loss
            p == TE && (α *= (ratio + (kcoa/n)^-2))
            γ += α
        end
        β = -im * γ
        if p == TE
           Z = k / β * η
        else
           Z = β / k * η
        end
        c.modes[q] = Mode(; n = mode.n, p, kcoa = mode.kcoa, γ, Z)
    end
    return c
end


function compute_kappa_matrix!(c1::CWG, c2::CWG)
    modes1 = c1.modes
    modes2 = c2.modes
    N1, N2 = length.((modes1, modes2))
    t = c2.a / c1.a  # Radius ratio
    if t < 1
        error("c1 radius greater than c2 radius")
    end

    kappas = zeros(Float64, N2, N1)  # Preallocation

    for q1 in 1:N1       # Loop over region I modes
        mode1 = modes1[q1]
        n1 = mode1.n   # radial mode index
        cn1 = mode1.kcoa # chi or chi' for Region 1
        for q2 in 1:N2     # Loop over region II modes
            mode2 = modes2[q2]
            n = mode2.n     # radial mode index
            cn2 = mode2.kcoa  # chi or chi' for Region 2
            if mode1.p == TE
                if mode2.p == TE
                    # Equation (19) of notes:
                    kappas[q2, q1] =
                        2 * cn2 * cn1 * Iint1(cn1, cn2 / t, J1chip[n1], J0chip[n1]) /
                        (t * J1chip[n] * J1chip[n1] * mysqrt(cn2^2 - 1) * mysqrt(cn1^2 - 1))
                elseif mode2.p == TM
                    # Equation (21) of reference:
                    kappas[q2, q1] =
                        2 * besselj1(cn2 / t) * J1chip[n1] /
                        (cn2 * mysqrt(cn1^2 - 1) * J2chi[n] * J1chip[n1])
                end
            elseif mode1.p == TM
                if mode2.p == TE
                    kappas[q2, q1] = 0.0  # Equation (20) of reference.
                #case TM
                elseif mode2.p == TM
                    # Equation (22) of reference:
                    kappas[q2, q1] =
                        2 * Iint1(cn1, cn2 / t, 0.0, J0chi[n1]) /
                        (t * J2chi[n1] * J2chi[n])
                end
            end
        end
    end

    return kappas

end


function Iint1(a::Float64, b::Float64, J1a::Float64, J0a::Float64)
    # the integral defined in Eq. (18) of the notes.
    if abs(a - b) < 1.0e-6
        integral = ((a^2 - 2) * (J1a)^2 + (a * J0a)^2) / (2 * a^2)
    else
        J1b = besselj1(b)
        integral =
            (a * J1a * (J1b / b - besselj0(b)) - b * J1b * (J1a / a - J0a)) /
            ((b - a) * (b + a))
    end
    return integral
end


function junction!(gsm::GSM, c1::CWG, c2::CWG, kappas::Matrix{Float64})
    # Compute scattering matrix entries for the junction of `c1` and `c2`. Ref. Eq. 15, 

    # Obtain number of modes on both sides of junction:
    Nm1 = length(c1.modes)
    Nm2 = length(c2.modes)

    # Use Eq. (15) to find P matrix:
    #P = Diagonal(rz2) * kappas * Diagonal(rz1)
    P = [complex(v) for v in kappas]
    for n in 1:Nm1
        rZ1 = sqrt(c1.modes[n].Z)
        for m in 1:Nm2
            P[m, n] *= rZ1
        end
    end
    for m in 1:Nm2
        rZ2inv = inv(sqrt(c2.modes[m].Z))
        for n in 1:Nm1
            P[m, n] *= rZ2inv
        end
    end

    # Find scattering matrix using Eq. (11) and the facts that 
    # 2I - (I - A) = I + A, and  2I + (A - I) = I + A
    PtP = transpose(P) * P
    #Srow1 = (I + PtP) \ [(I - PtP) 2*transpose(P)]
    RHS = [Iminus!(PtP) 2*transpose(P)]
    Srow1 = ldiv!(lu!(twoIminus!(PtP)), RHS)
    PPt = P * transpose(P)
    #Srow2 = (I + PPt) \ [2*P (PPt - I)]
    RHS = [2*P minusI!(PPt)]
    Srow2 = ldiv!(lu!(twoIplus!(PPt)), RHS)

    gsm.s11 .= @view Srow1[:, 1:Nm1]
    gsm.s12 .= @view Srow1[:, (Nm1 + 1):end]
    gsm.s21 .= @view Srow2[:, 1:Nm1]
    gsm.s22 .= @view Srow2[:, (Nm1 + 1):end]

    return gsm
end

function junction(c1::CWG, c2::CWG, kappas::Matrix{Float64})
    n1 = length(c1.modes)
    n2 = length(c2.modes)
    gsm = GSM(n1, n2)
    return junction!(gsm, c1::CWG, c2::CWG, kappas::Matrix{Float64})
end



"""
    propagate!(gsm, modes, l)

Modify a GSM to reflect propagation of the the modes through a length `l`.
## Arguments
- `gsm`: A `GSM` compatible with `γs` such that `size(gsm.s22) .== length(modes)`. On exit `gsm` is modified.
- `modes::AbstractVector{Mode}`: The modes to be propagated.  
- `l`: The length [inch] through which the modes should propagate/attenuate.
"""
function propagate!(a::GSM, modes::AbstractVector{Mode}, l::Real)
    n2 = size(a[2,2], 1)
    n2 ≠ length(modes) && @error "# modes not consistent" n2 length(modes) size(a[1,2]) exception = ErrorException
    #  Loop over each of the modes:
    for i in 1:n2
        p = exp(-modes[i].γ * l)
        a[1,2][:, i] .*= p
        a[2,1][i, :] .*= p
        a[2,2][:, i] .*= p
        a[2,2][i, :] .*= p
    end
    return a
end



"""
    cascade(c1::CWG, c2::CWG, kappas::Matrix=zeros(0,0)) -> (gsm::GSM, kappas)

Compute the GSM for the junction of two waveguide sections, including their lengths.

## Input Arguments
- `c1` and `c2`: The two circular waveguides.  It is assumed that the frequency-dependent 
  `modes` vector in each of these has been updated for the current analysis frequency.
- `kappas`: Optional kappa (frequency-independent mode coupling) matrix.  If empty (default),
  it will be computed and returned to the caller in the return value of the function.
"""
function cascade(c1::CWG, c2::CWG, kappas::Matrix=zeros(0,0))
    isempty(kappas) && (kappas = compute_kappa_matrix!(c1::CWG, c2::CWG))
    n1 = length(c1.modes)
    n2 = length(c2.modes)
    size(kappas) == (n2,n1) || error("Size error")
    gsm = GSM(n1, n1)
    propagate!(gsm, c1.modes, c1.l)
    gsm2 = junction(c1, c2, kappas)
    gsm = GSMs.cascade(gsm, gsm2)
    propagate!(gsm, c2.modes, c2.l)
    return gsm, kappas
end

end # module