using TE11CylCavity: CWG, setup_modes!, cascade
using BenchmarkTools

a1 = 0.2765
a2 = 0.3

for n in (2, 4, 8, 16, 32, 64, 128)
    cwg1 = CWG(a=a1, l=0.0)
    cwg2 = CWG(a=a2, l=0.0)

    n1 = n
    n2 = round(Int, cwg2.a/cwg1.a * n1)
    isodd(n2) && (n2 += 1)

    fghz = 21.0

    setup_modes!(cwg1, fghz, n1)
    setup_modes!(cwg2, fghz, n2)

    #kappas = Cav.CWGModeMatching.compute_kappa_matrix!(cwg1, cwg2)

    gsm, _ = cascade(cwg1,cwg2)

    #print("n1=$n1, n2=$n2, fghz=$fghz, ")
    s = gsm.s22[1,1]
    println("s22: ", round(abs(s), digits=5), ", ", round(rad2deg(angle(s)), digits=2))
    #println(round(abs(s), digits=4))
    #println(round(angle(s) |> rad2deg, digits=2))
end

# cascade take 1.5 ms for n1=54, n2=70