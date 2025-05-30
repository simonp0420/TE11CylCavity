# Here we adjust the cavity length and conductivity to agree with CST empty cavity result

using TE11CylCavity: TE11CylCavity as Cav



fghzstart = 15.3
σstart = 644_000 # Wall conductivity S/m
geom2 = (a = 0.35, d = 3.5, z1 = 0.5, z2 = 0.55)

fresgoal = 15.3894 # CST for empty cavity
Qgoal = 1801.3 # CST for empty cavity
(; σ, a, info) =  Cav.findσa(geom2, fresgoal, Qgoal, σstart)

println("σ and a optimized to produce same resonant frequency and Q as CST for empty cavity:")
@show σ
@show a

# Case 0: Empty Cavity
geom = merge(geom2, (; a))

ϵᵣ = 1.0
tanδ = 0.0
case0 = Cav.findfresQ(geom, σ, ϵᵣ, tanδ, fghzstart)

println()
println("With adjusted σ and a, S-parameter model prediction for empty cavity:")
@show case0

# Case 2:
#ϵᵣ = 5.7
#tanδ = 0.003
fresgoal = 15.34311 # CST for case2 dielectric sample
Qgoal = 1636.789 # CST for case2 dielectric sample
case2opt = Cav.findet(geom, σ, fresgoal, Qgoal)
println()
println("""Sample parameters predicted solely from CST resonant frequency and Q:
           (compare to actuals: ϵᵣ=5.7, tanδ=0.003)""")
@show case2opt


case2 = Cav.findfresQ(geom, σ, case2opt.ϵᵣ, case2opt.tanδ, fresgoal)
println()
println("Estimated dielectric properties produce these errors:")
@show (fresgoal, case2.fres - fresgoal)
@show (Qgoal, case2.Q - Qgoal)
nothing