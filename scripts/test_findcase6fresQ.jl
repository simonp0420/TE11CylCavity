# Case 6 from HFSS
println("Case 6 (HFSS Simulation of Kratos Cavity with Spacer)")
println("Actuals: ϵᵣ = 7.3, tanδ = 0.002\n")

# Define the geometry and provided measured or computed resonant frequency and Q for empty and with-sample cavity:
aI = 0.276 # Radius of cwg I (spacer) [inch]
aD = 0.3 # Radius of cwg D (dielectric sample) [inch]
aII = 0.3 # Radius of cwg II (final air region) [inch]
d = 2.99 # Total cavity length [inch]
lI = 0.2; lD = 0.05; lII = d - (lI + lD) # CWG lengths adding up to d = 3.5 inch
σstart = 789_889.0 # Inconel wall conductivity S/m
fQ0 = (; fres = 14.00063, Q = 1338.7645) # Provided values for empty cavity [GHz, unitless]
fQ = (; fres = 13.291171, Q = 1219.2253) # Provided values for cavity with dielectric sample

using MKL # This only for Intel CPUs
using TE11CylCavity: TE11CylCavity, CWG, setup_rect2cyl, setup_cyl2rect, setup_modes!, findfresQ, findσd, findet
 


# Define the rectangular-to-circular and circular-to-rectangular waveguide transitions:
setup_rect2cyl(pkgdir(TE11CylCavity, "s2pfiles", "Kratos_rwg2cwgII_transition.s2p"))
setup_cyl2rect(pkgdir(TE11CylCavity, "s2pfiles", "Kratos_cwgII2rwg_transition.s2p"))


dstart = d # Save for later use

cwgI = CWG(; a=aI, l=lI, σ=σstart)
cwgD = CWG(; a=aD, l=lD, σ=σstart)
cwgII = CWG(; a=aII, l=lII, σ=σstart)

if cwgI.a == cwgD.a == cwgII.a
    n1 = 1 # Only one mode required since no step in radius
else
    n1 = 120
end
n2 = ceil(Int, n1 * (cwgD.a / cwgI.a))

foreach((c, n) -> setup_modes!(c, fQ0.fres, n), (cwgI, cwgD, cwgII), (n1, n2, n2))

fresgoal = fQ0.fres # Provided for empty cavity
Qgoal = fQ0.Q # Provided for empty cavity

case0 = findfresQ(cwgI, cwgD, cwgII, fresgoal)
println("Initial run for empty cavity: ", case0)

# Here we adjust the cavity length and conductivity to agree with the provided empty cavity result:
cwgIadj, cwgDadj, cwgIIadj, _, _ =  findσd(cwgI, cwgD, cwgII, fresgoal, Qgoal)
σopt, dopt = cwgIadj.σ, sum(c -> c.l, (cwgIadj, cwgDadj, cwgIIadj))

println("σ and d optimized to produce same resonant frequency and Q as provided for empty cavity:")
println((;σstart, σopt))
println((; dstart, dopt))

caseadj = findfresQ(cwgIadj, cwgDadj, cwgIIadj, fresgoal)

println()
println("With adjusted σ and d, S-parameter model prediction for empty cavity:")
@show caseadj

etpredicted = TE11CylCavity.findet(cwgIadj, cwgDadj, cwgIIadj, fQ.fres, fQ.Q)
(; ϵᵣ, tanδ) = etpredicted
println()
println("Sample parameters predicted from provided with-sample resonant frequency and Q:")
println((; ϵᵣ, tanδ))

cwgDopt = CWG(; a = cwgDadj.a, l=cwgDadj.l, σ=cwgDadj.σ, modes=cwgDadj.modes, ϵᵣ=etpredicted.ϵᵣ, tanδ=etpredicted.tanδ)
casea = findfresQ(cwgIadj, cwgDopt, cwgIIadj, fresgoal)
if casea.fres > fresgoal + 0.4
    casea = findfresQ(cwgIadj, cwgDopt, cwgIIadj, fresgoal-0.2)
elseif casea.fres < fresgoal - 0.4
    casea = findfresQ(cwgIadj, cwgDopt, cwgIIadj, fresgoal+0.2)
end
println()
println("Estimated dielectric properties produce these errors:")
fresgoal = fQ.fres; Qgoal = fQ.Q
@show (fresgoal, casea.fres - fresgoal)
@show (Qgoal, casea.Q - Qgoal)


nothing