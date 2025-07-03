"""
    read_touchstone(filename) -> (; comments, FGHz, Smat)

Read the comments, frequencies, and S-parameter data from a Version 1 Touchstone file

## Input Arguments
- `filename::AbstractString`: The name of the Touchstone file. 

## Return Value: A named tuple with the following fields:
    - `comments::Vector{String}`: The comment lines extracted from the file
    - `FGHz::Vector{Float64}`: The frequencies in the file, converted (if necessary) to GHz.
    - `Smat::Array{ComplexF64, 3}`: The complex S-parameters.  `Smat` has dimensions `(np, np, nf)`
      where `np` is the number of ports and `nf = length(FGHz)` is the number of frequencies.
"""
function read_touchstone(filename::AbstractString)
    data = touchstone_load(filename)
    FGHz = data.f .* 1e-9
    comments = data.comments
    Smat = data.N
    return (; comments, FGHz, Smat)
end
    

"""
    findfresQ(filename) -> (fres, Q)

Find the resonant frequency and Q from a Touchstone file.

Assumes that the file contains 2-port S-parameters for a resonator and that the 
frequency sweep is wide enough to cover more than the 3 dB bandwidth of the device.
Uses cubic spline interpolation of the magnitude squared of S12 to find the relevant frequencies.

## Input Arguments
- `filename::AbstractString`: The name of the Touchstone file.

## Return Value
A named tuple with the following fields:
- `fres`: The interpolated resonant frequency.
- `Q`: The quality factor computed as the ratio of `fres` to the half-power bandwidth.  The half-power
  points on either side of `fres` are also found by cubic spline interpolation.
"""
function findfresQ(s2pfile::AbstractString)
    s2p = read_touchstone(s2pfile)
    FGHz = s2p.FGHz
    Smat = s2p.Smat
    np = size(Smat, 1)
    np == 2 || error("# ports = $np, must be 2")
    ptrans = abs2.(@view Smat[1,2,:])
    spl = Dierckx.Spline1D(FGHz, ptrans) # Cubic spline fit
    fmid = 0.5 * (first(FGHz) + last(FGHz))

    fres = Roots.find_zero((first(FGHz), last(FGHz)), Roots.Bisection()) do f
        return Dierckx.derivative(spl, f)
    end

    halfpower = 0.5 * spl(fres)

    f1 = Roots.find_zero((first(FGHz), fres), Roots.Bisection()) do f
        return spl(f) - halfpower
    end

    f2 = Roots.find_zero((fres, last(FGHz)), Roots.Bisection()) do f
        return spl(f) - halfpower
    end

    Q = fres / (f2 - f1)

    return (; fres, Q)
end