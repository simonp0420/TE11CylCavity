# Case2a: Here we adjust the cavity length and conductivity to agree with CST empty cavity result

using TE11CylCavity: TE11CylCavity as Cav

println("Case 2a")

fghzstart = 15.3
σstart = 644_000 # Wall conductivity S/m
geom2 = (a = 0.35, d = 3.5, z1 = 0.5, z2 = 0.55)

fresgoal = 15.3791 # CST for empty cavity
Qgoal = 1813.12 # CST for empty cavity
(; σ, d, info) =  Cav.findσd(geom2, fresgoal, Qgoal, σstart)

println("σ and d optimized to produce same resonant frequency and Q as CST for empty cavity:")

@show σ
@show d

# Case 0: Empty Cavity
geom = merge(geom2, (; d))

ϵᵣ = 1.0
tanδ = 0.0
case0 = Cav.findfresQ(geom, σ, ϵᵣ, tanδ, fresgoal)

println()
println("With adjusted σ and d, S-parameter model prediction for empty cavity:")
@show case0

# Case 2a:
#ϵᵣ = 5.7
#tanδ = 0.003
fresgoal = 15.3355  # CST for case2a dielectric sample
Qgoal = 1650.81 # CST for case2a dielectric sample
case2aopt = Cav.findet(geom, σ, fresgoal, Qgoal)
println()
println("""Sample parameters predicted solely from CST resonant frequency and Q:
           (compare to actuals: ϵᵣ=5.7, tanδ=0.003)""")
@show case2aopt


case2a = Cav.findfresQ(geom, σ, case2opt.ϵᵣ, case2opt.tanδ, fresgoal)
println()
println("Estimated dielectric properties produce these errors:")
@show (fresgoal, case2a.fres - fresgoal)
@show (Qgoal, case2a.Q - Qgoal)


nothing