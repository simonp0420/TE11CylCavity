# Case2a
println("Case 2a")
println("Actuals: ϵᵣ = 5.7, tanδ = 0.003\n")

# Define the geometry and provided measured or computed resonant frequency and Q for empty and with-sample cavity:
a = 0.35 # Common radius of all 3 sections of CWG [inch]
d = 3.5 # Total cavity length [inch]
lI = 0.5; lD = 0.05; lII = d - (lI + lD) # CWG lengths adding up to d = 3.5 inch
σstart = 644_000 # Inconel wall conductivity S/m
fQ0 = (; fres = 15.3791, Q = 1813.12) # Provided values for empty cavity
fQ = (; fres = 15.3355, Q = 1650.81) # Provided values for case2a with dielectric sample
# Define the rectangular-to-circular and circular-to-rectangular waveguide transitions:
setup_rect2cyl(pkgdir(TE11CylCavity, "s2pfiles", "case2a_rwgcwg_transition.s2p"))
setup_cyl2rect(pkgdir(TE11CylCavity, "s2pfiles", "case2a_cwgrwg_transition.s2p"))


using TE11CylCavity: TE11CylCavity, CWG, setup_rect2cyl, setup_cyl2rect, setup_modes!, findfresQ, findσd, findet
 
dstart = d # Save for later use

cwgI = CWG(; a, l=lI, σ=σstart)
cwgD = CWG(; a, l=lD, σ=σstart)
cwgII = CWG(; a, l=lII, σ=σstart)

if cwgI.a == cwgD.a == cwgII.a
    n1 = 1 # Only one mode required since no step in radius
else
    n1 = 20
end
n2 = ceil(Int, n1 * (cwgD.a / cwgI.a))

foreach((c, n) -> setup_modes!(c, fQ0.fres, n), (cwgI, cwgD, cwgII), (n1, n2, n2))

fresgoal = fQ0.fres # Provided for empty cavity
Qgoal = fQ0.Q # Provided for empty cavity

case0 = findfresQ(cwgI, cwgD, cwgII, fresgoal)
println("Initial run for empty cavity: ", case0)

# Here we adjust the cavity length and conductivity to agree with the provided empty cavity result:
cwgIadj, cwgDadj, cwgIIadj, _ =  findσd(cwgI, cwgD, cwgII, fresgoal, Qgoal)
σopt, dopt = cwgIadj.σ, sum(c -> c.l, (cwgIadj, cwgDadj, cwgIIadj))

println("σ and d optimized to produce same resonant frequency and Q as provided for empty cavity:")
println((;σstart, σopt))
println((; dstart, dopt))

case0adj = findfresQ(cwgIadj, cwgDadj, cwgIIadj, fresgoal)

println()
println("With adjusted σ and d, S-parameter model prediction for empty cavity:")
@show case0adj

etpredicted = findet(cwgIadj, cwgDadj, cwgIIadj, fQ.fres, fQ.Q)
println()
println("Sample parameters predicted from provided with-sample resonant frequency and Q:")
println(etpredicted)

cwgDopt = CWG(; a = cwgDadj.a, l=cwgDadj.l, σ=cwgDadj.σ, modes=cwgDadj.modes, ϵᵣ=etpredicted.ϵᵣ, tanδ=etpredicted.tanδ)
case2a = findfresQ(cwgIadj, cwgDopt, cwgIIadj, fresgoal)
println()
println("Estimated dielectric properties produce these errors:")
fresgoal = fQ.fres; Qgoal = fQ.Q
@show (fresgoal, case2a.fres - fresgoal)
@show (Qgoal, case2a.Q - Qgoal)


nothing