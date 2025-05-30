using TE11CylCavity: TE11CylCavity as Cav



fghzstart = 15.3
σ = 644_000 # Wall conductivity S/m
geom2 = (a = 0.35, d = 3.5, z1 = 0.5, z2 = 0.55)

# Case 0: Empty Cavity
ϵᵣ = 1.0
tanδ = 0.0
case0 = Cav.findfresQ(geom2, σ, ϵᵣ, tanδ, fghzstart)
@show case0

# Case 2:
ϵᵣ = 5.7
tanδ = 0.003
case2 = Cav.findfresQ(geom2, σ, ϵᵣ, tanδ, fghzstart)
@show case2

nothing