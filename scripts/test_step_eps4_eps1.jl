using TE11CylCavity: CWG, setup_modes!, cascade
using BenchmarkTools

a1 = 0.276
a2 = 0.3

for n in (2, 4, 8, 16, 32, 64, 128)
    cwg1 = CWG(a = a1, l = 0.0)
    cwg2 = CWG(a = a2, ϵᵣ = 4.0, l = 0.05 / 2)
    cwg3 = CWG(a = a2, l = 0.0)

    n1 = n
    n2 = round(Int, cwg2.a/cwg1.a * n1)
    isodd(n2) && (n2 += 1)

    fghz = 18.0

    for (cwg, n) in ((cwg1, n1), (cwg2, n2), (cwg3, n2))
        setup_modes!(cwg, fghz, n)
    end

    gsm1, _ = cascade(cwg1,cwg2)
    gsm2, _ = cascade(cwg2, cwg3)
    gsm = cascade(gsm1, gsm2)
    #print("n1=$n1, n2=$n2, fghz=$fghz, ")
    s = gsm.s21[1,1]
    #println("s11: ", round(abs(s), digits=5), ", ", round(rad2deg(angle(s)), digits=2))
    println(round(abs(s), digits=4))
    #println(round(angle(s) |> rad2deg, digits=2))
end
