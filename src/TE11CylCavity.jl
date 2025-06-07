module TE11CylCavity

using StaticArrays: @SMatrix
using PRIMA: newuoa
using Dierckx: Dierckx # Spline1D and derivative
using Roots: Roots  # find_zero

include("CWGModeMatching.jl")
using .CWGModeMatching: CWGModeMatching as cmm

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
function cascade(a, b)
    den = 1 - a[2,2] * b[1,1]  
    c = @SMatrix [a[1,1]+a[1,2]*a[2,1]*b[1,1]/den   a[1,2]*b[2,1]/den
                  a[2,1]*b[1,2]/den                 b[2,2]+b[1,2]*b[2,1]*a[2,2]/den]
    return c
end


let s11r, s11i, s12r, s12i, s22r, s22i, flow, fhigh, fname # Establish a new hard scope 

    global srect2cyl, setup_transition

    """
        setup_transition(filename::AbstractString)

    Read in the 2-port scattering parameters from the Touchstone file specified 
    by `filename`.  Set up for cubic spline interpolation of the scattering parameters
    over frequency to be performed by `srect2cyl`

    """
    function setup_transition(filename::AbstractString)
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
        srect2cyl(f)
    Return the 2x2 S matrix of the lossy transition from rectangular (port 1) to circular 
    (port 2) waveguide at frequency `f` (in GHz) obtained via cubic spline interpolation into
    the S-parameter data established by by `setup_transition`.
    """
    function srect2cyl(f)
        flow ≤ f ≤ fhigh || 
          @warn """Requested frequency $f is outside interval [$flow,$fhigh] covered by $fname)
                   Continuing with extrapolation. This warning will only be displayed once"""  maxlog=1
        s11 = complex(s11r(f), s11i(f))
        s12 = complex(s12r(f), s12i(f))
        s22 = complex(s22r(f), s22i(f))
        S = @SMatrix [s11 s12; s12 s22]  
        #S = @SMatrix [0.998358*cis(deg2rad(177.2))   0.015691*cis(deg2rad(87.14))
        #               0.015691*cis(deg2rad(87.14))   0.998563*cis(deg2rad(178.48))]
        return S
    end

end # let


                       
"""
    sim_cavity_s21(geom, Rs, ϵᵣ, tanδ, fghz) -> smat

Simulate the value of S₁₂ for a cylindrical cavity with dielectric sample inserted.

## Input Arguments
- `geom`: A named tuple with the following fields:
  * `d`: The cavity length (between lossy shorting plates) in inches.
  * `a`: The cavity radius in inches.
  * `z1`, `z2`: The locations of the inserted dielectric sample boundaries, where `0 ≤ z1 < z2 ≤ d`.
- `Rs`: The cavity wall surface resistance in ohms.
- `ϵᵣ, tanδ`: The dielectric constant and loss tangent of the inserted dielectric sample.
- `fghz`: The frequency in GHz.

## Return Value
- `smat`: The complex scattering matrix measured between the two rectangular waveguides at each end.
"""                       
function sim_cavity_smat(geom, Rs, ϵᵣ, tanδ, fghz)
    cz = zero(ComplexF64)
    rootepsr = sqrt(ϵᵣ)

    # Region 1 (vacuum) parameters:
    k1 = kcalc(fghz)           # wavenumber (radians/inch)
    η1 = η₀
    β1 = βcalc(k1, geom.a, Rs/η1)  # Guide propagation constant (radians/inch)
    #Zh1 = k1 / β1 * η1        # Modal impedance

    # Region 2 (dielectric sample) parameters:
    k2 = kcalc(fghz, ϵᵣ, tanδ)
    η2 = η1 / rootepsr
    β2 = βcalc(k2, geom.a, real(Rs/η2))
    #Zh2 = k2 / β2 * η2  # Modal impedance
    
    # Begin cascade calculations...
    Sr2c = srect2cyl(fghz) # Rectangular to circular wg transition
    S = Sr2c

    # First Region 1 interval:
    t = cis(-β1 * geom.z1)
    b = @SMatrix [cz t; t cz]
    S = cascade(S, b)

    # Region 1->2 Interface Plane:
    den = β1 + β2
    b11 = (β1 - β2) / den
    b12 = 2 * sqrt(β1) * sqrt(β2) / den # Two sqrts to avoid subtle bug in complex domain
    b = @SMatrix [b11 b12; b12 -b11]
    S = cascade(S, b)

    # Region 2 propagation:
    t = cis(-β2 * (geom.z2 - geom.z1))
    b = @SMatrix [cz t; t cz]
    S = cascade(S, b)

    # Region 2->1 Interface Plane (note sign changes of s11, s22):
    b = @SMatrix [-b11 b12; b12 b11]
    S = cascade(S, b)

    # Final Region 1 interval:
    t = cis(-β1 * (geom.d - geom.z2))
    b = @SMatrix [cz t; t cz]
    S = cascade(S, b)

    # Transition from circular to rectangular waveguide:
    b = @SMatrix [Sr2c[2,2] Sr2c[2,1]; Sr2c[1,2] Sr2c[1,1]]
    S = cascade(S, b)
    return S
end

"""
    findfres(geom, σ, ϵᵣ, tanδ, fstart; rhoend=1e-12) -> fres, info

Find the cavity resonant frequency where |S₁₂| is a maximum.

## Positional Arguments
- `geom`: A named tuple with the following fields:
  * `d`: The cavity length (between lossy shorting plates) in inches.
  * `a`: The cavity radius in inches.
  * `z1`, `z2`: The locations of the inserted dielectric sample boundaries, where `0 ≤ z1 < z2 ≤ d`.
- `σ`: The bulk conductivity of the cavity wall in MKS units [S/m].
- `ϵᵣ, tanδ`: The dielectric constant and loss tangent of the inserted dielectric sample.
- `fstart`: The initial guess for the resonant frequency in GHz.

## Optional Keyword Arguments
- `rhoend`: Used as the `rhoend` argument to the `PRIMA.newuoa` function. 

## Return Values
- `fres`: An estimate for the resonant frequency in GHz.
- `info`: The `info` struct returned by `PRIMA.newuoa`.
"""
function findfres(geom, σ, ϵᵣ, tanδ, fstart; rhoend=1e-12)
    Rsstart = σ2Rs(σ, fstart)
    Sstart = sim_cavity_smat(geom, Rsstart, ϵᵣ, tanδ, fstart)
    tstart = abs(Sstart[1,2])
    x, info = newuoa([0.0]; rhoend) do x
        fghz = 1e-3 * x[1] * fstart + fstart # Unscale the objective variable
        Rs = σ2Rs(σ, fghz)
        S = sim_cavity_smat(geom, Rs, ϵᵣ, tanδ, fghz)
        s12abs = abs(S[1,2])
        return tstart / s12abs
    end
    fopt = 1e-3 * x[1] * fstart + fstart
    return fopt, info
end

"""
    findf1(geom, σ, ϵᵣ, tanδ, fres; rhoend=1e-12) -> f1, info

Find the cavity lower half-power frequency where |S₁₂|² is 1/2 of its maximum (resonance) value.

## Positional Arguments
- `geom`: A named tuple with the following fields:
  * `d`: The cavity length (between lossy shorting plates) in inches.
  * `a`: The cavity radius in inches.
  * `z1`, `z2`: The locations of the inserted dielectric sample boundaries, where `0 ≤ z1 < z2 ≤ d`.
- `σ`: The bulk conductivity of the cavity wall in MKS units [S/m].
- `ϵᵣ, tanδ`: The dielectric constant and loss tangent of the inserted dielectric sample.
- `fres`: The resonant frequency in GHz, where |S₁₂| takes its maximum value.

## Optional Keyword Arguments
- `rhoend`: Used as the `rhoend` argument to the `PRIMA.newuoa` function. 

## Return Values
- `f1`: An estimate for the lower 1/2 power frequency in GHz.
- `info`: The `info` struct returned by `PRIMA.newuoa`.
"""
function findf1(geom, σ, ϵᵣ, tanδ, fres; rhoend=1e-12)
    Rsres = σ2Rs(σ, fres)
    Sres = sim_cavity_smat(geom, Rsres, ϵᵣ, tanδ, fres)
    tgoal = inv(sqrt(2)) * abs(Sres[1,2])
    x, info = newuoa([0.0]; rhoend) do x
        fghz = (1 - 1e-6 * x[1]^2) * fres
        Rs = σ2Rs(σ, fghz)
        S = sim_cavity_smat(geom, Rs, ϵᵣ, tanδ, fghz)
        s12abs = abs(S[1,2])
        return (s12abs - tgoal)^2
    end
    fopt = (1 - 1e-6 * x[1]^2) * fres
    return fopt, info
end

"""
    findf2(geom, σ, ϵᵣ, tanδ, fres; rhoend=1e-12) -> f2, info

Find the cavity upper half-power frequency where |S₁₂|² is 1/2 of its maximum (resonance) value.

## Positional Arguments
- `geom`: A named tuple with the following fields:
  * `d`: The cavity length (between lossy shorting plates) in inches.
  * `a`: The cavity radius in inches.
  * `z1`, `z2`: The locations of the inserted dielectric sample boundaries, where `0 ≤ z1 < z2 ≤ d`.
- `σ`: The bulk conductivity of the cavity wall in MKS units [S/m].
- `ϵᵣ, tanδ`: The dielectric constant and loss tangent of the inserted dielectric sample.
- `fres`: The resonant frequency in GHz, where |S₁₂| takes its maximum value.

## Optional Keyword Arguments
- `rhoend`: Used as the `rhoend` argument to the `PRIMA.newuoa` function. 

## Return Values
- `f2`: An estimate for the upper 1/2 power frequency in GHz.
- `info`: The `info` struct returned by `PRIMA.newuoa`.
"""
function findf2(geom2, σ, ϵᵣ, tanδ, fres; rhoend=1e-12)
    Rsres = σ2Rs(σ, fres)
    Sres = sim_cavity_smat(geom2, Rsres, ϵᵣ, tanδ, fres)
    tgoal = inv(sqrt(2)) * abs(Sres[1,2])
    x, info = newuoa([0.0]; rhoend) do x
        fghz = (1 + 1e-6 * x[1]^2) * fres
        Rs = σ2Rs(σ, fghz)
        S = sim_cavity_smat(geom2, Rs, ϵᵣ, tanδ, fghz)
        s12abs = abs(S[1,2])
        return (s12abs - tgoal)^2
    end
    fopt = (1 + 1e-6 * x[1]^2) * fres
    return fopt, info
end


"""
    findfresQ(geom, σ, ϵᵣ, tanδ, fstart) -> (; fres, Q)

Find the resonant frequency and Q of the given model.

## Input Arguments
- `geom`: A named tuple with the following fields:
  * `d`: The cavity length (between lossy shorting plates) in inches.
  * `a`: The cavity radius in inches.
  * `z1`, `z2`: The locations of the inserted dielectric sample boundaries, where `0 ≤ z1 < z2 ≤ d`.
- `σ`: The cavity wall conductivity in Siemen/meter.
- `ϵᵣ, tanδ`: The dielectric constant and loss tangent of the inserted dielectric sample.
- `fstart`: The initial guess at resonant frequency in GHz.
"""
function findfresQ(geom, σ, ϵᵣ, tanδ, fstart)
    fres, infores = findfres(geom, σ, ϵᵣ, tanδ, fstart)
    f1, info1 = findf1(geom, σ, ϵᵣ, tanδ, fres)
    f2, info2 = findf2(geom, σ, ϵᵣ, tanδ, fres)
    Q = fres / (f2 - f1)
    return (; fres, Q)
end

"""
findσd(geom, fresgoal, Qgoal, σstart) -> (; σ, d, info)

Adjust conductivity and cylinder length (slightly) to yield the desired resonant frequency and Q of empty cavity.

The idea is to slightly adjust the wall conductivity and cavity length so that the S-parameter model produces
exactly the same resonant frequency and Q as measured.  These updated values can then be used when seeking the
dielectric constant and loss tangent corresponding to measured resonant frequency and Q with the sample inserted.

## Input Arguments
- `geom`: A named tuple with the following fields:
  * `d`: The cavity length (between lossy shorting plates) in inches.
  * `a`: The cavity radius in inches.
  * `z1`, `z2`: The locations of the inserted dielectric sample boundaries, where `0 ≤ z1 < z2 ≤ d`.
    For this function these are immaterial, since the analysis will assume vacuum electrical parameters in the sample region.
- `fresgoal`: The desired empty cavity resonant frequency in GHz.
- `Qgoal`: The desired empty cavity quality factor.
- `σstart`: An initial guess for the cavity wall conductivity in Siemen/meter.

## Return Value
A named tuple with the following fields:
- `σ`: The optimized conductivity [S/m].
- `d`: The optimized cavity length [inch].
- `info`: The `info` struct returned by `PRIMA.newuo`
"""
function findσd(geom, fresgoal, Qgoal, σstart)
    rhobeg = 0.1
    rhoend = 1e-12
    dstart = geom.d
    x, info = newuoa([1.0, 1.0]; rhobeg, rhoend) do x
        ϵᵣ = 1.0
        tanδ = 0.0
        σ = x[1]^2 * σstart
        d = x[2]^2 * dstart
        geomtest = merge(geom, (; d))
        (; fres, Q) = findfresQ(geomtest, σ, ϵᵣ, tanδ, fresgoal)
        objective = sqrt(1e6 * (fres - fresgoal)^2 + 1e-6 * (Q - Qgoal)^2)
        return objective
    end
    σ = x[1]^2 * σstart
    d = x[2]^2 * dstart
    return (; σ, d, info)

end

"""
findσa(geom, fresgoal, Qgoal, σstart) -> (; σ, a, info)

Adjust conductivity and cylinder radius (slightly) to yield the desired resonant frequency and Q of empty cavity.

The idea is to slightly adjust the wall conductivity and cavity radius so that the S-parameter model produces
exactly the same resonant frequency and Q as measured.  These updated values can then be used when seeking the
dielectric constant and loss tangent corresponding to measured resonant frequency and Q with the sample inserted.

## Input Arguments
- `geom`: A named tuple with the following fields:
  * `d`: The cavity length (between lossy shorting plates) in inches.
  * `a`: The cavity radius in inches.
  * `z1`, `z2`: The locations of the inserted dielectric sample boundaries, where `0 ≤ z1 < z2 ≤ d`.
    For this function these are immaterial, since the analysis will assume vacuum electrical parameters in the sample region.
- `fresgoal`: The desired empty cavity resonant frequency in GHz.
- `Qgoal`: The desired empty cavity quality factor.
- `σstart`: An initial guess for the cavity wall conductivity in Siemen/meter.

## Return Value
A named tuple with the following fields:
- `σ`: The optimized conductivity [S/m].
- `a`: The optimized cavity radius [inch].
- `info`: The `info` struct returned by `PRIMA.newuo`
"""
function findσa(geom, fresgoal, Qgoal, σstart)
    rhobeg = 0.1
    rhoend = 1e-12
    astart = geom.a
    x, info = newuoa([1.0, 1.0]; rhobeg, rhoend) do x
        ϵᵣ = 1.0
        tanδ = 0.0
        σ = x[1]^2 * σstart
        a = x[2]^2 * astart
        geomtest = merge(geom, (; a))
        (; fres, Q) = findfresQ(geomtest, σ, ϵᵣ, tanδ, fresgoal)
        objective = sqrt(1e6 * (fres - fresgoal)^2 + 1e-6 * (Q - Qgoal)^2)
        return objective
    end
    σ = x[1]^2 * σstart
    a = x[2]^2 * astart
    return (; σ, a, info)

end


"""
    findet(geom, σ, fresgoal, Qgoal) -> (; ϵᵣ, tanδ, info)

Determine sample dielectric constan and loss tangent to yield the desired resonant frequency and Q cavity with sample.

## Input Arguments
- `geom`: A named tuple with the following fields:
  * `d`: The cavity length (between lossy shorting plates) in inches.
  * `a`: The cavity radius in inches.
  * `z1`, `z2`: The locations of the inserted dielectric sample boundaries, where `0 ≤ z1 < z2 ≤ d`.
- `σ`: The wall conductivity in S/m.
- `fresgoal`: The desired empty cavity resonant frequency in GHz.
- `Qgoal`: The desired empty cavity quality factor.

## Return Value
A named tuple with the following fields:
- `ϵᵣ`: An estimate of the dielectric constant of the sample.
- `tanδ`: An estimate of the loss tangent of the sample.
- `info`: The `info` struct returned by `PRIMA.newuo`
"""
function findet(geom, σ, fresgoal, Qgoal)
    rhobeg = 1.0
    rhoend = 1e-13
    x, info = newuoa([0.0, 0.0]; rhobeg, rhoend) do x
        ϵᵣ = 1.0 + x[1]^2
        tanδ = x[2]^2
        (; fres, Q) = findfresQ(geom, σ, ϵᵣ, tanδ, fresgoal)
        objective = sqrt(1e6 * (fres - fresgoal)^2 + 1e-6 * (Q - Qgoal)^2)
        return objective
    end
    ϵᵣ = 1.0 + x[1]^2
    tanδ = x[2]^2
    return (; ϵᵣ, tanδ, info)
end

end # module