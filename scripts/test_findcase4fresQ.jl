# Case4: Kratos Cavity Case Without Spacer
#Here we adjust the cavity length and conductivity to agree with CST empty cavity result
geom = (a = 0.3, d = 2.99, z1 = 0.2, z2 = 0.25)
fQ0 = (; fres = 16.50816, Q = 1661.18) # CST for empty cavity
fQ = (; fres = 15.89933, Q = 1659.412) # CST for case4 dielectric sample
σstart = 789889 # Inconel wall conductivity S/m
rwg2cavii_s2pfile = joinpath(@__DIR__, "KratosCavityRWG2CWGII.s2p") # RWG/CWG transition S-parameters

using TE11CylCavity: TE11CylCavity as Cav

Cav.setup_transition(rwg2cavii_s2pfile)

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

σopt = σstart
dopt = geom.d

# Case 0: Empty Cavity
geomopt = merge(geom, (; d=dopt))

case0adj = Cav.findfresQ(geomopt, σopt, 1.0, 0.0, fresgoal)

println()
println("With adjusted σ and d, S-parameter model prediction for empty cavity:")
@show case0adj


etpredicted = Cav.findet(geomopt, σopt, fQ.fres, fQ.Q)
etpredicted = (; ϵᵣ = round(etpredicted.ϵᵣ, digits=3), tanδ = round(etpredicted.tanδ, sigdigits=3))
println()
println("Sample parameters predicted solely from CST resonant frequency and Q:")
println(etpredround)


casenew = Cav.findfresQ(geomopt, σopt, etpredicted.ϵᵣ, etpredicted.tanδ, fresgoal)
println()
println("Estimated dielectric properties produce these errors:")
@show (fQ.fres, fQ.fres - casenew.fres)
@show (fQ.Q, fQ.Q - casenew.Q)


nothing