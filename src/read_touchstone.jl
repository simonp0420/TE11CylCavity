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
    comments = String[]
    FGHz = Float64[]
    ports = 0 # Initialize count of ports
    ports_from_ext = false
    row = 0   # Initialize row counter
    line = ""
    kline = 0
    GHz_multiplier = 1.0
    tocomplex = (a, b) -> a * cis(deg2rad(b))
    
    lines = readlines(filename)

    # Compute # ports from file extension:
    ext = lowercase(splitext(filename)[2])
    if length(ext) > 3 && ext[1:2] == ".s" && ext[end] == 'p'&& all(isdigit, ext[3:end-1])
        ports = parse(Int, ext[3:end-1])
        ports_from_ext = true
    end
        
    
    kline = 0
    while kline < length(lines)
        kline += 1
        line = lines[kline]
        # Read comments
        if line[1] == '!'
            push!(comments,line)
            if !ports_from_ext
                k = findfirst("Port[", line)
                if !isnothing(k)
                    k2 = findnext(']', line, k[end]+1)
                    ports = max(ports, parse(Int, line[k[end]+1:k2[begin]-1]))
                end
            end
        elseif line[1] == '#'  # Format line
            line = uppercase(line)
            larray = split(line)
            if larray[2] == "GHZ"
                GHz_multiplier = 1.0
            elseif larray[2] == "HZ"
                GHz_multiplier = 1.0e-9
            elseif larray[2] == "MHZ"
                GHz_multiplier = 1.0e-3
            else
                error("Unknown frequency units specifier: $(larray[2])")
            end
            if !isequal(larray[3], "S")
                error("Third entry of the format line was \"$(larray[3])\" rather than \"S\" ")
            end
            if larray[4] == "MA" 
                tocomplex = (a,b) -> a * cis(b*pi/180)
            elseif larray[4] == "DB"
                tocomplex = (a,b) -> 10^(a/20) * cis(b*pi/180)
            elseif larray[4] == "RI"
                tocomplex = (a,b) -> complex(a,b)
            else
                error("Unknown format specifier \"$(larray[4])\" ")
            end
        else 
            break  # Data line found
        end
    end
    # Done with header comments.  lines[kline] contains first row of S-parameter data and the first frequency
    # Find the empty lines:
    empty = findall(isempty, lines)
    if !isempty(empty)
        Nf = length(empty)  # Number of frequencies
    else
        Nf = length(lines) - kline + 1
    end
    
    Smat = zeros(ComplexF64, (ports, ports, Nf))
    kline -= 1
    for kf = 1:Nf
        dat = Float64[]
        while length(dat) < 2*ports*ports + 1
            kline += 1
            append!(dat, parse.(Float64, split(lines[kline])))
        end
        push!(FGHz, GHz_multiplier * dat[1])
        k = 0
        for col = 1:ports
            for row = 1:ports
                k += 1
                Smat[row,col,kf] = tocomplex(dat[2k],dat[2k+1])
            end
        end
        
        # Skip any comments following this block:
        while kline < length(lines)
            kline += 1
            line = lines[kline]
            if isempty(line)
                continue
            elseif line[1] == '!' || line[1] == '#'
                push!(comments,line)
                continue
            else
                kline -= 1
                break
            end
        end
    end
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