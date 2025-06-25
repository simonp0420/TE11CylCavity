# Case2a: Here we adjust the cavity length and conductivity to agree with CST empty cavity result

using TE11CylCavity: CWG, setup_rect2cyl, setup_cyl2rect, setup_modes!, findfresQ, findσd, findet
 
a = 0.35 # Common radius of all 3 sections of CWG
d = 3.5 # Total length
dstart = d # Save for later use
lI = 0.5; lD = 0.05; lII = d - (lI + lD) # CWG lengths adding up to d = 3.5 inch
σstart = 644_000 # Inconel wall conductivity S/m

cwgI = CWG(; a, l=lI, σ=σstart)
cwgD = CWG(; a, l=lD, σ=σstart)
cwgII = CWG(; a, l=lII, σ=σstart)

n1 = 20
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

case0 = findfresQ(cwgI, cwgD, cwgII, fresgoal)
println("Initial run for empty cavity: ", case0)

t =  findσd(cwgI, cwgD, cwgII, fresgoal, Qgoal)
cwgIadj = t.cwgI
cwgDadj = t.cwgD
cwgIIadj = t.cwgII
σopt, dopt = cwgIadj.σ, sum(c -> c.l, (cwgIadj, cwgDadj, cwgIIadj))

println("σ and d optimized to produce same resonant frequency and Q as CST for empty cavity:")
println((;σstart, σopt))
println((; dstart, dopt))

case0adj = findfresQ(cwgIadj, cwgDadj, cwgIIadj, fresgoal)

println()
println("With adjusted σ and d, S-parameter model prediction for empty cavity:")
@show case0adj

# Case 2a:
#ϵᵣ = 5.7
#tanδ = 0.003
etpredicted = findet(cwgIadj, cwgDadj, cwgIIadj, fQ.fres, fQ.Q)
println()
println("""Sample parameters predicted solely from CST resonant frequency and Q:
           (compare to actuals: ϵᵣ=5.7, tanδ=0.003)""")
println(etpredicted)

cwgDopt = CWG(; a = cwgDadj.a, l=cwgDadj.l, σ=cwgDadj.σ, modes=cwgDadj.modes, ϵᵣ=etpredicted.ϵᵣ, tanδ=etpredicted.tanδ)
case2a = findfresQ(cwgIadj, cwgDopt, cwgIIadj, fresgoal)
println()
println("Estimated dielectric properties produce these errors:")
fresgoal = fQ.fres; Qgoal = fQ.Q
@show (fresgoal, case2a.fres - fresgoal)
@show (Qgoal, case2a.Q - Qgoal)


nothing