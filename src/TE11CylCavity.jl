module TE11CylCavity

export CWG, setup_modes!, cascade, junction, junction!, propagate!, setup_rect2cyl, setup_cyl2rect, scyl2rect, srect2cyl

using StaticArrays: @SMatrix
using PRIMA: newuoa, bobyqa
using Dierckx: Dierckx # Spline1D and derivative
using Roots: Roots  # find_zero
using Touchstone: touchstone_load

include("CWGModeMatching.jl")
import .CWGModeMatching: cascade  # To be extended
using .CWGModeMatching: CWG, setup_modes!, junction, junction!, propagate!, compute_kappa_matrix

include("Constants.jl")
using .Constants: c₀mks, c₀, μ₀, ϵ₀, η₀

include("read_touchstone.jl")


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
    kcalc(f, ϵᵣ = 1.0, tanδ = 0.0)

Compute the complex intrinsic wavenumber at frequency `f` (in GHz) of a medium with dielectric constant `ϵᵣ`
and loss tangent `tanδ`.  The returned value will be in units of radians/inch.
"""
function kcalc(f, ϵᵣ = 1.0, tanδ = 0.0)
    λ₀ = c₀ / f
    k0 = 2π / λ₀
    k = k0 * sqrt(ϵᵣ * complex(1.0, -tanδ))
    imag(k) > 0 && (k = -k)
    return k
end


"""
    kcocalc(a)
Compute cutoff radial wavenumber for the TE11 mode in a circular waveguide of radius `a`.
Units will be inverse length in whatever length units `a` was expressed in.
"""
function kcocalc(a)
    p′₁₁ = 1.84118378134065930264363 # First zero of J1'
    return p′₁₁ / a
end

"""
    βcalc(k, a, Rsnorm) -> β
Compute the complex propagation constant β for the TE11 mode in a circular waveguide of radius `a`.
## Input Arguments
- `k`: The medium instrinsic wavenumber (radians/inch). 
- `a`: The radius of the cylindrical waveguide (inches).
- `Rsnorm`: The ratio Rs/η, where Rs is the metal surface resistance and η is the intrinsic impedance of the medium.
## Return Value
- `β`: Complex propagation constant [radians/inch] such that `cis(-β*z)` is the z-dependence of a positive-going wave.
"""
function βcalc(k, a, Rsnorm)
    kco = kcocalc(a)
    βdiel = sqrt(k^2 - kco^2)
    imag(βdiel) > 0 && (βdiel = -βdiel)
    kco = kcocalc(a)
    p′11 = kcocalc(1.0)
    kratio = (kco / k)^2 |> real
    α = Rsnorm / a * inv(sqrt(1 - kratio)) * (kratio + inv(p′11^2 - 1)) # Collin FTOGW Table 5.5
    β = βdiel + complex(0.0, -α)
    return β
end


"""
    cascade(a, b) -> c
Cascade two-port circuits represented by the 2x2 scattering matrices `a` and `b`, and return result in `c`.
Assumes that port 2 of `a` is connected to port 1 of `b`. 
"""
function cascade(a::AbstractMatrix, b::AbstractMatrix)
    size(a) == size(b) == (2,2) || throw(ArgumentError("a and b must be 2×2 matrices"))
    den = 1 - a[2,2] * b[1,1]  
    c = @SMatrix [a[1,1]+a[1,2]*a[2,1]*b[1,1]/den   a[1,2]*b[2,1]/den
                  a[2,1]*b[1,2]/den                 b[2,2]+b[1,2]*b[2,1]*a[2,2]/den]
    return c
end


let s11r, s11i, s12r, s12i, s22r, s22i, flow, fhigh, fname # Establish a new hard scope 

    global srect2cyl, setup_rect2cyl

    """
        setup_rect2cyl(filename::AbstractString)

    Read in the 2-port scattering parameters from the Touchstone file specified 
    by `filename`.  Set up for cubic spline interpolation of the scattering parameters
    over frequency to be performed by `srect2cyl`

    """
    function setup_rect2cyl(filename::AbstractString)
        (; comments, FGHz, Smat) = read_touchstone(filename)
        flow, fhigh = first(FGHz), last(FGHz)
        s11r = Dierckx.Spline1D(FGHz, real.(@view Smat[1,1,:]), k=3, bc="extrapolate")
        s12r = Dierckx.Spline1D(FGHz, real.(@view Smat[1,2,:]), k=3, bc="extrapolate")
        s22r = Dierckx.Spline1D(FGHz, real.(@view Smat[2,2,:]), k=3, bc="extrapolate")
        s11i = Dierckx.Spline1D(FGHz, imag.(@view Smat[1,1,:]), k=3, bc="extrapolate")
        s12i = Dierckx.Spline1D(FGHz, imag.(@view Smat[1,2,:]), k=3, bc="extrapolate")
        s22i = Dierckx.Spline1D(FGHz, imag.(@view Smat[2,2,:]), k=3, bc="extrapolate")
    end

    """
        srect2cyl(f)
    Return the 2x2 S matrix of the lossy transition from rectangular (port 1) to circular 
    (port 2) waveguide at frequency `f` (in GHz) obtained via cubic spline interpolation into
    the S-parameter data established by by `setup_rect2cyl`.  It is assumed that this is
    the end of the cavity where the spacer is inserted, so the circular waveguide radius
    is that of the spacer inner radius.
    """
    function srect2cyl(f)
        flow ≤ f ≤ fhigh || 
          @warn """Requested frequency $f is outside interval [$flow,$fhigh] covered by $fname)
                   Continuing with extrapolation. This warning will only be displayed once"""  maxlog=1
        s11 = complex(s11r(f), s11i(f))
        s12 = complex(s12r(f), s12i(f))
        s22 = complex(s22r(f), s22i(f))
        S = @SMatrix [s11 s12; s12 s22]  
        return S
    end

end # let



let s11r, s11i, s12r, s12i, s22r, s22i, flow, fhigh, fname # Establish a new hard scope 

    global scyl2rect, setup_cyl2rect

    """
        setup_cyl2rect(filename::AbstractString)

    Read in the 2-port scattering parameters from the Touchstone file specified 
    by `filename`.  Set up for cubic spline interpolation of the scattering parameters
    over frequency to be performed by `scyl2rect`.  

    """
    function setup_cyl2rect(filename::AbstractString)
        fname = filename
        (; comments, FGHz, Smat) = read_touchstone(filename)
        flow, fhigh = first(FGHz), last(FGHz)
        s11r = Dierckx.Spline1D(FGHz, real.(@view Smat[1,1,:]), k=3, bc="extrapolate")
        s12r = Dierckx.Spline1D(FGHz, real.(@view Smat[1,2,:]), k=3, bc="extrapolate")
        s22r = Dierckx.Spline1D(FGHz, real.(@view Smat[2,2,:]), k=3, bc="extrapolate")
        s11i = Dierckx.Spline1D(FGHz, imag.(@view Smat[1,1,:]), k=3, bc="extrapolate")
        s12i = Dierckx.Spline1D(FGHz, imag.(@view Smat[1,2,:]), k=3, bc="extrapolate")
        s22i = Dierckx.Spline1D(FGHz, imag.(@view Smat[2,2,:]), k=3, bc="extrapolate")
    end

    """
        scyl2rect(f)
    Return the 2x2 S matrix of the lossy transition from circular (port 1) to rectangular
    (port 2) waveguide at frequency `f` (in GHz) obtained via cubic spline interpolation into
    the S-parameter data established by by `setup_cyl2rect`.   It is assumed that this is
    the end of the cavity without the spacer, so the circular waveguide radius
    is that of the larger cavity.
    """
    function scyl2rect(f)
        flow ≤ f ≤ fhigh || 
          @warn """Requested frequency $f is outside interval [$flow,$fhigh] covered by $fname)
                   Continuing with extrapolation. This warning will only be displayed once"""  maxlog=1
        s11 = complex(s11r(f), s11i(f))
        s12 = complex(s12r(f), s12i(f))
        s22 = complex(s22r(f), s22i(f))
        S = @SMatrix [s11 s12; s12 s22]  
        return S
    end

end # let


                       
"""
    sim_cavity_smat!(cwgI, cwgD, cwgII, kappasID, kappasDII, fghz) -> smat

Simulate the value of S₁₂ for a cylindrical cavity with dielectric sample inserted.

Assumes that `setup_rect2cyl` and `setup_cyl2rect` have already been called to initialize data needed by 
`scyl2rect` and `srect2cyl`.

## Input Arguments
- `cwgI::CWG`: A `CWG` describing the Region I circular waveguide (i.e., the spacer).
- `cwgD::CWG`: A `CWG` describing the dielectric sample region of circular waveguide.
- `cwgII::CWG`: A `CWG` describing the Region II circular waveguide (i.e., the cavity minus the spacer and dielectric region).
- `kappas`:  Frequency-independent mode coupling matrix for the junction of `cwgI` and `cwgD`.
- `fghz`: The frequency in GHz.

It is assumed that `cwgI`, `cwgII`, and `cwgD` have their `modes` field initialized, although they are not assumed to be 
updated for the current frequency specified in `fghz`.  These fields will then be modified.

    ## Return Value
- `smat`: The complex 2×2 scattering matrix measured between the two rectangular waveguides at each end of the cavity.
"""                       
function sim_cavity_smat!(cwgI::CWG, cwgD::CWG, cwgII::CWG, kappasID::Matrix{Float64}, kappasDII::Matrix{Float64}, fghz)
    foreach(cwg -> setup_modes!(cwg, fghz), (cwgI, cwgD, cwgII)) # Update the modes for this frequency
    # Multimode cascades:
    SID, _ = cascade(cwgI, cwgD, kappasID)
    SDII = junction(cwgD, cwgII, kappasDII)
    propagate!(SDII, cwgII.modes, cwgII.l)
    S = cascade(SID, SDII)
    # Single-mode cascades:
    S2x2 = @SMatrix [S[1,1][1,1] S[1,2][1,1] ; S[2,1][1,1] S[2,2][1,1]] # Retain only TE10 mode
    Sr2c = srect2cyl(fghz) # Rectangular to circular wg transition
    S2x2 = cascade(Sr2c, S2x2)
    Sc2r = scyl2rect(fghz) # Circular to rectangular wg transition
    S2x2 = cascade(S2x2, Sc2r)
    return S2x2
end

"""
    findfres(cwgI::CWG, cwgD::CWG, cwgII::CWG, fstart; kappasID, kappasDII, rhoend=1e-12) -> fres, info

Find the cavity resonant frequency where |S₁₂| is a maximum.

## Positional Arguments
- `cwgI`, `cwgD`, `cwgII`: `CWG` instances describing the Region I, Dielectric Sample, and Region II waveguides. The `modes`
  field of each should be initialized.
- `fstart`: The initial guess for the resonant frequency in GHz.

## Optional Keyword Arguments
- `kappasID`, `kappasDII`: Frequency-independent coupling matrices for the Region I/Dielectric Sample and Dielectric Sample/Region II
  interfaces, respectively.  If not passed they will be computed, at a cost of execution time.
- `rhoend`: Used as the `rhoend` argument to the `PRIMA.newuoa` function. 

## Return Values
- `fres`: An estimate for the resonant frequency in GHz.
- `info`: The `info` struct returned by `PRIMA.newuoa`.
"""
function findfres(
    cwgI::CWG,
    cwgD::CWG,
    cwgII::CWG,
    fstart;
    kappasID = compute_kappa_matrix(cwgI, cwgD),
    kappasDII = compute_kappa_matrix(cwgD, cwgII),
    rhoend = 1e-12)
    
    Sstart = sim_cavity_smat!(cwgI, cwgD, cwgII, kappasID, kappasDII, fstart)
    tstart = abs(Sstart[1,2])
    x, info = bobyqa([0.0]; rhoend) do x
        fghz = 1e-3 * x[1] * fstart + fstart # Unscale the objective variable
        S = sim_cavity_smat!(cwgI, cwgD, cwgII, kappasID, kappasDII, fghz)
        s12abs = abs(S[1,2])
        return tstart / s12abs
    end
    fopt = 1e-3 * x[1] * fstart + fstart
    return fopt, info
end

"""
    findf1(cwgI, cwgD, cwgII, fres; kappasID, kappasDII, rhoend=1e-12) -> f1, info

Find the cavity lower half-power frequency where |S₁₂|² is 1/2 of its maximum (resonance) value.

## Positional Arguments
- `cwgI`, `cwgD`, `cwgII`: `CWG` instances describing the Region I, Dielectric Sample, and Region II waveguides. The `modes`
  field of each should be initialized.
- `fres`: The resonant frequency in GHz, where |S₁₂| takes its maximum value.

## Optional Keyword Arguments
- `kappasID`, `kappasDII`: Frequency-independent coupling matrices for the Region I/Dielectric Sample and Dielectric Sample/Region II
  interfaces, respectively.  If not passed they will be computed, at a cost of execution time.
- `rhoend`: Used as the `rhoend` argument to the `PRIMA.newuoa` function. 

## Return Values
- `f1`: An estimate for the lower 1/2 power frequency in GHz.
- `info`: The `info` struct returned by `PRIMA.newuoa`.
"""
function findf1(
    cwgI::CWG,
    cwgD::CWG,
    cwgII::CWG,
    fres;
    kappasID = compute_kappa_matrix(cwgI, cwgD),
    kappasDII = compute_kappa_matrix(cwgD, cwgII),
    rhoend=1e-10)

    Sres = sim_cavity_smat!(cwgI, cwgD, cwgII, kappasID, kappasDII, fres)
    tgoal = inv(sqrt(2)) * abs(Sres[1,2])

    Qmax = 2e4
    Qmin = 2e2
    rhobeg = inv(2Qmax)
    xl = [1 - inv(2Qmin)]
    xu = [1 - inv(2Qmax)]
    x, info = bobyqa([xu[1] - rhobeg]; xl, xu, rhobeg, rhoend) do x
        fghz = x[1] * fres
        S = sim_cavity_smat!(cwgI, cwgD, cwgII, kappasID, kappasDII, fghz)
        s12abs = abs(S[1,2])
        return (s12abs - tgoal)^2
    end
    fopt = x[1] * fres
    return fopt, info
end

"""
    findf2(cwgI, cwgD, cwgII, fres; kappasID, kappasDII, rhoend=1e-12) -> f2, info

Find the cavity upper half-power frequency where |S₁₂|² is 1/2 of its maximum (resonance) value.

## Positional Arguments
- `cwgI`, `cwgD`, `cwgII`: `CWG` instances describing the Region I, Dielectric Sample, and Region II waveguides. The `modes`
  field of each should be initialized.
- `fres`: The resonant frequency in GHz, where |S₁₂| takes its maximum value.

## Optional Keyword Arguments
- `kappasID`, `kappasDII`: Frequency-independent coupling matrices for the Region I/Dielectric Sample
  and Dielectric Sample/Region II interfaces, respectively.  If not passed they will be computed, costing execution time.
- `rhoend`: Used as the `rhoend` argument to the `PRIMA.newuoa` function. 

## Return Values
- `f2`: An estimate for the upper 1/2 power frequency in GHz.
- `info`: The `info` struct returned by `PRIMA.newuoa`.
"""
function findf2(
    cwgI::CWG,
    cwgD::CWG,
    cwgII::CWG,
    fres;
    kappasID = compute_kappa_matrix(cwgI, cwgD),
    kappasDII = compute_kappa_matrix(cwgD, cwgII),
    rhoend=1e-10)

    Sres = sim_cavity_smat!(cwgI, cwgD, cwgII, kappasID, kappasDII, fres)
    tgoal = inv(sqrt(2)) * abs(Sres[1,2])

    Qmax = 2e4
    Qmin = 2e2
    rhobeg = inv(2Qmax)
    xl = [1 + inv(2Qmax)]
    xu = [1 + inv(2Qmin)]
    x, info = bobyqa([xl[1] + rhobeg]; xl, xu, rhobeg, rhoend) do x
        fghz = x[1] * fres
        S = sim_cavity_smat!(cwgI, cwgD, cwgII, kappasID, kappasDII, fghz)
        s12abs = abs(S[1,2])
        return (s12abs - tgoal)^2
    end
    fopt = x[1] * fres
    return fopt, info
end


"""
    findfresQ(cwgI, cwgD, cwgII, fstart; kappasID, kappasDII) -> (; fres, Q)

Find the resonant frequency and Q of the given model.

## Positional Input Arguments
- `cwgI`, `cwgD`, `cwgII`: `CWG` instances describing the Region I, Dielectric Sample, and Region II waveguides. The `modes`
  field of each should be initialized.
- `fstart`: The initial guess at resonant frequency in GHz.
## Optional Keyword Arguments
- `kappasID`, `kappasDII`: Frequency-independent coupling matrices for the Region I/Dielectric Sample
  and Dielectric Sample/Region II interfaces, respectively.  If not passed they will be computed, costing execution time.
"""
function findfresQ(
    cwgI,
    cwgD,
    cwgII,
    fstart;
    kappasID = compute_kappa_matrix(cwgI, cwgD),
    kappasDII = compute_kappa_matrix(cwgD, cwgII)
    )
    
    fres, infores = findfres(cwgI, cwgD, cwgII, fstart; kappasID, kappasDII)
    f1, info1 = findf1(cwgI, cwgD, cwgII, fres; kappasID, kappasDII)
    f2, info2 = findf2(cwgI, cwgD, cwgII, fres; kappasID, kappasDII)
    Q = fres / (f2 - f1)
    return (; fres, Q)
end

"""
findσd(cwgI, cwgD, cwgII, fresgoal, Qgoal) -> (; cwgI, cwgD, cwgII, info)

Adjust conductivity of all circular waveguides and cwgII cylinder length (slightly) to yield the desired resonant frequency and Q of empty cavity.

The idea is to slightly adjust the wall conductivity and cavity length so that the S-parameter model produces
exactly the same resonant frequency and Q as measured.  These updated values can then be used when seeking the
dielectric constant and loss tangent corresponding to measured resonant frequency and Q with the sample inserted.

## Input Arguments
- `cwgI`, `cwgD`, `cwgII`: `CWG` instances describing the Region I, Dielectric Sample, and Region II waveguides. Each of these
  should specify the same conductivity (`σ` field) value. Also the `modes` field of each should be initialized.
  `cwgD` should be initialized for a dielectric constant of unity and zero loss tangent.
- `fresgoal`: The desired empty cavity resonant frequency in GHz.
- `Qgoal`: The desired empty cavity quality factor.

## Return Values
A named tuple with the following fields:
- `cwgI`, `cwgD`, `cwgII`: Same as the input values except that the conductivity (σ) values will be updated with the new
  value and the length (`l` field) of `cwgII` will be adjusted as well.
- `d`: The optimized cavity length [inch].
- `infolII`: The `info` struct returned by `PRIMA.newuo` when optimizing `cwgII.l`
- `infoσ`: The `info` struct returned by `PRIMA.newuo` when optimizing σ for all 3 guides.
"""
function findσd(cwgI, cwgD, cwgII, fresgoal, Qgoal)
    kappasID = compute_kappa_matrix(cwgI, cwgD)
    kappasDII = compute_kappa_matrix(cwgD, cwgII)
    σstart = cwgI.σ
    lIIstart = cwgII.l

    # For d
    rhobeg = 0.1
    rhoend = 1e-8
    x, infolII = bobyqa([1.0]; rhobeg, rhoend, xl=[0.8], xu=[1.2]) do x
        lII = x[1] * lIIstart
        cwgIItest = CWG(; a=cwgII.a, l=lII, σ = cwgII.σ, modes=cwgII.modes)
        fres, _ = findfres(cwgI, cwgD, cwgIItest, fresgoal; kappasID, kappasDII)
        objective = (fres/fresgoal - 1)^2
        return objective
    end
    lII = x[1] * lIIstart

    # For σ
    rhobeg = 0.1
    rhoend = 1e-8
    x, infoσ = bobyqa([1.0]; rhobeg, rhoend, xl=[0.8], xu=[1.2]) do x
        σ = x[1] * σstart
        cwgItest = CWG(; a=cwgI.a, l=cwgI.l, σ, modes=cwgI.modes)
        cwgDtest = CWG(; a=cwgD.a, l=cwgD.l, σ, modes=cwgD.modes)
        cwgIItest = CWG(; a=cwgII.a, l=lII, σ, modes=cwgII.modes)
        f1, _ = findf1(cwgItest, cwgDtest, cwgIItest, fresgoal; kappasID, kappasDII)
        f2, _ = findf2(cwgItest, cwgDtest, cwgIItest, fresgoal; kappasID, kappasDII)
        Q = fresgoal / (f2 - f1)
        objective = (Q/Qgoal - 1)^2
        return objective
    end
    σ = x[1] * σstart
    cwgI = CWG(; a=cwgI.a, l=cwgI.l, σ, modes=cwgI.modes)
    cwgD = CWG(; a=cwgD.a, ϵᵣ=cwgD.ϵᵣ, tanδ=cwgD.tanδ, l=cwgD.l, σ, modes=cwgD.modes)
    cwgII = CWG(; a=cwgII.a, l=lII, σ, modes=cwgII.modes)

    return (; cwgI, cwgD, cwgII, infolII, infoσ)
end



"""
    findet(cwgI, cwgD, cwgII, fresgoal, Qgoal; n1=120) -> (; ϵᵣ, tanδ, info)

Determine sample dielectric constant and loss tangent to yield the desired resonant frequency and Q cavity with sample.

## Positional Input Arguments
- `cwgI`, `cwgD`, `cwgII`: `CWG` instances describing the Region I, Dielectric Sample, and Region II waveguides. Each of these
  should specify the same conductivity (`σ` field) value. The `modes` field of these should be initialized to an empty array.
- `fresgoal`: The desired empty cavity resonant frequency in GHz.
- `Qgoal`: The desired empty cavity quality factor.
## Optional Keyword Arguments
- `n1`: The number of modes used in `cwgI` for the multi-mode GSM analysis.  Defaults to 150.

## Return Value
A named tuple with the following fields:
- `ϵᵣ`: An estimate of the dielectric constant of the sample.
- `tanδ`: An estimate of the loss tangent of the sample.
- `infoϵᵣ`: The `info` struct returned by `PRIMA.bobyqa` when optimizing `cwgD.ϵᵣ`.
- `infotanδ`: The `info` struct returned by `PRIMA.bobyqa` when optimizing `cwgD.tanδ`.
"""
function findet(cwgI::CWG, cwgD::CWG, cwgII::CWG, fresgoal, Qgoal; n1=120)

    n1 = iszero(length(cwgI.modes)) ? n1 : length(cwgI.modes)
    setup_modes!(cwgI, fresgoal, n1)
    nD = ceil(Int, n1 * cwgD.a / cwgI.a)
    nD = iszero(length(cwgD.modes)) ? nD : length(cwgD.modes)
    setup_modes!(cwgD, fresgoal, nD)
    n2 = ceil(Int, n1 * cwgII.a / cwgI.a)
    n2 = iszero(length(cwgII.modes)) ? n2 : length(cwgII.modes)
    setup_modes!(cwgII, fresgoal, n2)
    
    kappasID = compute_kappa_matrix(cwgI, cwgD)
    kappasDII = compute_kappa_matrix(cwgD, cwgII)

    # Starting values:
    ϵᵣ = 2.5
    tanδ = 10^(-2.5)

    # For dielectric constant
    rhobeg = 0.35
    rhoend = 1e-6
    x, infoϵᵣ = bobyqa([ϵᵣ]; rhobeg, rhoend, xl=[1.0], xu=[15.0]) do x
        ϵᵣ = x[1]
        cwgDtest = CWG(; a=cwgD.a, l=cwgD.l, ϵᵣ, tanδ, σ=cwgD.σ, modes=cwgD.modes)
        fres, _ = findfres(cwgI, cwgDtest, cwgII, fresgoal; kappasID, kappasDII)
        objective = (fres/fresgoal - 1)^2
        return objective
    end
    ϵᵣ = x[1]
    
    # For loss tangent
    rhobeg = 0.2
    rhoend = 1e-8
    x, infotanδ = bobyqa([-log10(tanδ)]; rhobeg, rhoend, xl=[0.5], xu=[5.0]) do x
        tanδ = 10^(-x[1])
        cwgDtest = CWG(; a=cwgD.a, l=cwgD.l, ϵᵣ, tanδ, σ=cwgD.σ, modes=cwgD.modes)
        f1, info1 = findf1(cwgI, cwgDtest, cwgII, fresgoal; kappasID, kappasDII)
        f2, info2 = findf2(cwgI, cwgDtest, cwgII, fresgoal; kappasID, kappasDII)
        Q = fresgoal / (f2 - f1)
        objective = (Q/Qgoal - 1)^2
        return objective
    end
    tanδ = 10^(-x[1])

    return (; ϵᵣ, tanδ, infoϵᵣ, infotanδ)
end

end # module