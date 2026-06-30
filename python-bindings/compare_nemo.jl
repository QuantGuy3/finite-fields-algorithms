# compare_nemo.jl — Nemo.jl (FLINT) reference for the alcp showcase.
#
# Reads identical inputs to those alcp uses (so the comparison is apples-to-
# apples), multiplies and factors over GF(p), and writes results + timings.
#   julia compare_nemo.jl <in.txt> <out.txt>
#
# in.txt:   line1 p_mult | line2 A coeffs | line3 B coeffs | line4 p_factor | line5 F coeffs
# out.txt:  MULT_MS, MULT, FACTOR_MS, FACTOR (canonical monic factor multiset)
using Nemo

parse_csv(s) = isempty(strip(s)) ? Int[] : parse.(Int, split(strip(s), ","))

# integer coefficients (little-endian) of a GF(p) polynomial
coeffs_le(h) = degree(h) < 0 ? Int[] : [Int(lift(ZZ, coeff(h, i))) for i in 0:degree(h)]

function main()
    inp, outp = ARGS[1], ARGS[2]
    L = readlines(inp)
    pmul = parse(Int, strip(L[1]))
    A = parse_csv(L[2]);  B = parse_csv(L[3])
    pfac = parse(Int, strip(L[4]))
    Fc = parse_csv(L[5])

    out = IOBuffer()
    reps = 5

    # ---- multiplication over GF(pmul) ----
    Rm, _ = polynomial_ring(Nemo.Native.GF(pmul), "x")
    fa = Rm(A); fb = Rm(B)
    prod = fa * fb                      # warm-up (also JIT-compiles)
    t = @elapsed for _ in 1:reps; prod = fa * fb end
    println(out, "MULT_MS ", t * 1000 / reps)
    println(out, "MULT ", join(coeffs_le(prod), ","))

    # ---- factorization over GF(pfac) ----
    Rf, _ = polynomial_ring(Nemo.Native.GF(pfac), "x")
    ff = Rf(Fc)
    fac = factor(ff)                    # warm-up
    t2 = @elapsed for _ in 1:reps; fac = factor(ff) end
    items = String[]
    for (g, e) in fac
        push!(items, join(coeffs_le(g), ",") * ":" * string(e))
    end
    println(out, "FACTOR_MS ", t2 * 1000 / reps)
    println(out, "FACTOR ", join(items, ";"))

    write(outp, String(take!(out)))
end

main()
