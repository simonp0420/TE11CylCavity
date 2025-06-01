# Here we adjust the cavity length and conductivity to agree with CST empty cavity result

using TE11CylCavity: TE11CylCavity as Cav



fghzstart = 15.3
σstart = 644_000 # Wall conductivity S/m
geom2 = (a = 0.35, d = 3.5, z1 = 0.5, z2 = 0.55)

fresgoal = 15.3894 # CST for empty cavity
Qgoal = 1801.3 # CST for empty cavity
(; σ, d, info) =  Cav.findσd(geom2, fresgoal, Qgoal, σstart)

println("σ and d optimized to produce same resonant frequency and Q as CST for empty cavity:")
@show σ
@show d

# Case 0: Empty Cavity
geom = merge(geom2, (; d))

ϵᵣ = 1.0
tanδ = 0.0
case0 = Cav.findfresQ(geom, σ, ϵᵣ, tanδ, fghzstart)

println()
println("With adjusted σ and d, S-parameter model prediction for empty cavity:")
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

# Examine sensitivity of ϵᵣ and tanδ to changes in the resonant frequency and Q of the second cavity
using FiniteDifferences: grad, central_fdm
fun_eps_f(f) = Cav.findet(geom, σ, f, Qgoal).ϵᵣ
fun_eps_Q(Q) = Cav.findet(geom, σ, fresgoal, Q).ϵᵣ
∂ϵᵣ∂f = central_fdm(11, 1)(fun_eps_f, fresgoal)[1] # -88.73020968895703
∂ϵᵣ∂Q = central_fdm(11, 1; factor=2e18)(fun_eps_Q, Qgoal)[1] # -3.475112954414834e-6

fun_tand_f(f) = Cav.findet(geom, σ, f, Qgoal).tanδ
fun_tand_Q(Q) = Cav.findet(geom, σ, fresgoal, Q).tanδ
∂tanδ∂f = central_fdm(3, 1; factor=1e12)(fun_tand_f, fresgoal)[1] # 0.1457377324694493
∂tanδ∂Q = central_fdm(11, 1; factor=1e17)(fun_tand_Q, Qgoal)[1] # -4.3673063289868904e-5

s = """[∂ϵᵣ/∂f      ∂ϵᵣ/∂Q  
         ∂tanδ/∂f   ∂tanδ/∂Q ]"""
println(s)
nothing