# Here we adjust the frequency and Q goals based on the ratio of the S-param model to CST values for empty cavity
# This produces crummy results:
#Sample parameters predicted solely from CST resonant frequency and Q:
#(compare to actuals: ϵᵣ=5.7, tanδ=0.003)
#case2opt = (ϵᵣ = 6.423452439760667, tanδ = 0.0020796954511283976, info = PRIMA.Info(2.278946403967382e-5, 185, PRIMA.SMALL_TR_RADIUS, 0.0, Float64[], Float64[]))

using TE11CylCavity: TE11CylCavity as Cav



fghzstart = 15.3
σstart = 644_000 # Wall conductivity S/m
geom2 = (a = 0.35, d = 3.5, z1 = 0.5, z2 = 0.55)

ϵᵣ = 1.0
tanδ = 0.0
case0 = Cav.findfresQ(geom2, σ, ϵᵣ, tanδ, fghzstart)

# CST results for empty cavity
fresCSTcase0 = 15.3894
QCSTcase0 = 1801.3

fresCSTcase2 = 15.34311
QCSTcase2 = 1636.789 # CST for case2 dielectric sample

# Optimization goals:
fresgoal = fresCSTcase2 * (case0.fres / fresCSTcase0)
Qgoal = QCSTcase2 * (case0.Q / QCSTcase0)

# Case 2:
#ϵᵣ = 5.7
#tanδ = 0.003
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