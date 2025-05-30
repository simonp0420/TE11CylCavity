using TE11CylCavity: TE11CylCavity as Cav



fghz = 15.38
λ0 = Cav.c₀mks / (fghz * 1e9) # meters
k0 = 2π / λ0 # radians/meter
σ = 644_000 # S/m
Rs = sqrt(k0 * Cav.η₀ / 2σ)
geom2 = (a = 0.35, d = 3.5, z1 = 0.5, z2 = 0.55)

# Case 0: Empty Cavity
ϵᵣ = 1.0
tanδ = 0.0
case0 = Cav.findfresQ(geom2, Rs, ϵᵣ, tanδ, fghz)
@show case0

# Case 2:
ϵᵣ = 5.7
tanδ = 0.003
case2 = Cav.findfresQ(geom2, Rs, ϵᵣ, tanδ, fghz)
@show case2

nothing