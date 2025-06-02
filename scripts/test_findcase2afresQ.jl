# Case2a: Here we adjust the cavity length and conductivity to agree with CST empty cavity result
geom = (a = 0.35, d = 3.5, z1 = 0.5, z2 = 0.55)
fQ0 = (; fres = 15.3791, Q = 1813.12) # CST for empty cavity
fQ = (; fres = 15.3355, Q = 1650.81) # CST for case2a dielectric sample
σstart = 644_000 # Inconel wall conductivity S/m

using TE11CylCavity: TE11CylCavity as Cav

Cav.setup_transition(joinpath(@__DIR__, "rwgcwg_transition.s2p"))
println("Case 2a")

fghzstart = fQ0.fres

fresgoal = fQ0.fres # CST for empty cavity
Qgoal = fQ0.Q # CST for empty cavity
case0 = Cav.findfresQ(geom, σstart, 1.0, 0.0, fresgoal)
println("Initial run for empty cavity: ", case0)

t =  Cav.findσd(geom, fresgoal, Qgoal, σstart)
σopt, dopt = t.σ, t.d

println("σ and d optimized to produce same resonant frequency and Q as CST for empty cavity:")
dstart = geom.d
println((;σstart, σopt))
println((; dstart, dopt))

# Case 0: Empty Cavity
geomopt = merge(geom, (; d=dopt))

case0adj = Cav.findfresQ(geomopt, σopt, 1.0, 0.0, fresgoal)

println()
println("With adjusted σ and d, S-parameter model prediction for empty cavity:")
@show case0adj

# Case 2a:
#ϵᵣ = 5.7
#tanδ = 0.003
etpredicted = Cav.findet(geomopt, σopt, fQ.fres, fQ.Q)
println()
println("""Sample parameters predicted solely from CST resonant frequency and Q:
           (compare to actuals: ϵᵣ=5.7, tanδ=0.003)""")
println(etpredicted)


case2a = Cav.findfresQ(geomopt, σopt, etpredicted.ϵᵣ, etpredicted.tanδ, fresgoal)
println()
println("Estimated dielectric properties produce these errors:")
@show (fresgoal, case2a.fres - fresgoal)
@show (Qgoal, case2a.Q - Qgoal)


nothing