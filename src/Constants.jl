module Constants

export c₀mks, c₀, μ₀, ϵ₀, η₀

"Speed of light [m/s]"
const c₀mks = 299792458.0

"Speed of light [inch*GHz]"
const c₀ = c₀mks / 25.4e6 

"Permeability of free space [H/m]"
const μ₀ = 1.25663706212 * 1e-6

"Permeability of free space [F/m]"
const ϵ₀ = 1 / (μ₀ * c₀mks^2)

"Intrinsic impedance of free space [F/m]"
const η₀ = sqrt(μ₀ / ϵ₀)

end # module