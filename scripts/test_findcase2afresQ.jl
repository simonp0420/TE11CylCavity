# Case2a: Here we adjust the cavity length and conductivity to agree with CST empty cavity result

using TE11CylCavity: CWG, setup_rect2cyl, setup_cyl2rect, setup_modes!, findfresQ, findσd, findet
 
a = 0.35 # Common radius of all 3 sections of CWG
d = 3.5 # Total length
lI = 0.2; lD = 0.05; lII = d - (lI + lD) # CWG lengths adding up to d = 3.5 inch
σstart = 644_000 # Inconel wall conductivity S/m

cwgI = CWG(; a, l=lI, σ=σstart)
cwgD = CWG(; a, l=lD, σ=σstart)
cwgII = CWG(; a, l=lII, σ=σstart)

n1 = 2
n2 = ceil(Int, n1 * cwgD.a / cwgI.a); isodd(n2) && (n2 += 1)

fQ0 = (; fres = 15.3791, Q = 1813.12) # CST for empty cavity
fQ = (; fres = 15.3355, Q = 1650.81) # CST for case2a dielectric sample

foreach((c, n) -> setup_modes!(c, fQ0.fres, n), (cwgI, cwgD, cwgII), (n1, n2, n2))



setup_rect2cyl(joinpath(@__DIR__, "case2a_rwgcwg_transition.s2p"))
setup_cyl2rect(joinpath(@__DIR__, "case2a_cwgrwg_transition.s2p"))
println("Case 2a")

fghzstart = fQ0.fres

fresgoal = fQ0.fres # CST for empty cavity
Qgoal = fQ0.Q # CST for empty cavity

using TE11CylCavity: sim_cavity_smat!, compute_kappa_matrix
kappasID = compute_kappa_matrix(cwgI, cwgD)
kappasDII = compute_kappa_matrix(cwgD, cwgII)
@show sim_cavity_smat!(cwgI, cwgD, cwgII, kappasID, kappasDII, fQ0.fres)

case0 = findfresQ(cwgI, cwgD, cwgII, fresgoal)
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