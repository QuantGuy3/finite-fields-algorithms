"""
PRÁCTICA ALCP 2

Grupo: Ismael Amador, Ángela Ruiz, Miquel Sandonís, Elia Torres y Claudia Vicario
"""


# para esta actividad necesitamos la anterior (cuerpos_finitos.py)
# completa, pero quitando las funciones de factorización en anillo_fp_x
# y en anillo_fq_x, ya que las implementaremos aquí (no se olviden de
# quitarlas)
import cuerpos_finitos as cf
from copy import deepcopy
import random

# square-free factorization
# input: fpx --> anillo_fp_x
# input: f --> polinomio fabricado por fpx (objeto opaco) no nulo
# output: g = producto de los factores irreducibles mónicos de f, es decir,
# si f = c * f1^e1 * f2^e2 * ... * fr^er con los fi irreducibles mónicos
# distintos entre si, ei >= 1, c en fp, entonces g = f1 * f2 * ... * fr
# square-free factorization para anillo_fp_x (característica p)
def sqfree_fact_fpx(fpx, f):
    F = fpx._to_list(f)
    p = fpx.p

    # derivada
    if len(F) <= 1:
        df = [0]
    else:
        df = [0] * (len(F)-1)
        for i in range(1, len(F)):
            df[i-1] = (i * F[i]) % p
    df = fpx._trim(df)

    # caso derivada nula: f(x) = g(x^p). Construir g y recursar.
    if fpx._gradoL(df) == -1:
        g_coef = [ F[i] for i in range(0, len(F), p) ]   # coeficientes de g en orden creciente
        g_poly = fpx._from_list(fpx._trim(g_coef))
        return sqfree_fact_fpx(fpx, g_poly)

    # caso derivada no nula: gcd y tomar la parte sin multiplicidades
    cL = fpx._gcdL(F, df)
    qL, _ = fpx._divmodL(F, cL)   # qL = f / gcd(f,f')
    gL = fpx._trim(qL)

    # limpiar g para que sea mónico
    if fpx._gradoL(gL) > 0:
        lead = gL[fpx._gradoL(gL)]
        if lead != 1:
            inv_lead = fpx.K._inv(lead)
            gL = [fpx.K._mul(coef, inv_lead) for coef in gL]

    # eliminar f de todos los factores encontrados en g
    residuoL = list(F)
    while True:
        yL = fpx._gcdL(residuoL, gL)
        if fpx._gradoL(yL) <= 0:
            break
        residuoL, _ = fpx._divmodL(residuoL, yL)

    # lo que queda en residuoL es una potencia p-ésima
    if fpx._gradoL(residuoL) > 0:
        h_coef = [0] * ((len(residuoL) - 1) // p + 1)
        for i in range(0, len(residuoL), p):
            h_coef[i // p] = residuoL[i]

        h_poly = fpx._from_list(fpx._trim(h_coef))

        g2_poly = sqfree_fact_fpx(fpx, h_poly)

        g_poly = fpx._from_list(gL)

                # Unimos los factores de ambas partes evitando duplicados (LCM)
        producto = g_poly * g2_poly
        interseccion = fpx.gcd(g_poly, g2_poly)
        # producto / interseccion = LCM(g_poly, g2_poly)
        return fpx.div(producto, interseccion)

    else:

        return fpx._from_list(gL)


# distinct-degree factorization
# input: fpx --> anillo_fp_x
# input: g --> polinomio de fpx (objeto opaco) que es producto de factores
# irreducibles mónicos distintos cada uno con multiplicidad uno
# output: [h1, h2, ..., hr], donde hi = producto de los factores irreducibles
# mónicos de h de grado = i, el último hr debe ser no nulo y por supuesto
# g = h1 * h2 * ... * hr
def didegr_fact_fpx(fpx, g):
    p = fpx.p
    xL = [0, 1]
    gL = fpx._to_list(g)

    H_list = []
    H_current = fpx._divmodL(xL, gL)[1]
    i = 1

    deg_g = fpx._gradoL(gL)

    while i<= deg_g and deg_g>0:
        H_current = fpx._pot_modL(H_current, p, gL)
        dL = fpx._gcdL(gL, fpx._sumaL(H_current, fpx._negL(xL)))

        if not fpx._is_unoL(dL):
            H_list.append(fpx._from_list(dL))
            qL, _ = fpx._divmodL(gL, dL)
            gL = qL
            deg_g = fpx._gradoL(gL)
            H_current = fpx._divmodL(H_current, gL)[1] if deg_g >= 0 else [0]

        else:
            H_list.append(fpx.uno())

        i += 1

        if i > deg_g // 2 and deg_g > 0:
            while len(H_list) < deg_g -1:
                H_list.append(fpx.uno())
            H_list.append(fpx._from_list(gL))
            break

    return H_list

# equal-degree factorization
# input: fpx --> anillo_fp_x (supondremos p impar)
# input: r --> int
# input: h --> polinomio de fpx (objeto opaco) que es producto de factores
# irreducibles mónicos distintos de grado r con multiplicidad uno
# output: [u1, ..., us], donde h = u1 * u2* ... * us y los ui son irreducibles
# mónicos de grado = r
def eqdegr_fact_fpx(fpx, r, h):
    p = fpx.p
    hL = fpx._to_list(h)
    deg_h = fpx._gradoL(hL)

    if deg_h <= r:
        if deg_h > 0:
            return [h]
        else:
            return []
    
    factors = []
    stack = [hL]

    while stack: 
        gL = stack.pop()
        deg_g = fpx._gradoL(gL)

        if deg_g == r:
            if deg_g > 0:
                factors.append(fpx._from_list(gL))
            continue

        if deg_g <= 0:
            continue
        
        found_split = False
        for _ in range(max(10, deg_g)):
            coefs = [random.randrange(p) for _ in range(deg_g)]
            coefs.append(0)
            aL = fpx._trim(coefs)
            if fpx._gradoL(aL) <=0:
                continue

            exp = (pow(p, r) -1)//2
            bL = fpx._pot_modL(aL, exp, gL)

            dL = fpx._gcdL(gL, fpx._sumaL(bL, fpx._negL([1])))

            dr = fpx._gradoL(dL)

            if 0 < dr < deg_g:
                qL, _ = fpx._divmodL(gL, dL)
                stack.append(dL)
                stack.append(qL)
                found_split = True
                break

        if not found_split:
            if fpx._gradoL(gL) > 0:
                factors.append(fpx._from_list(gL))  

    return factors

# multiplicidad de factor irreducible mónico
# input: fpx --> anillo_fp_x
# input: f --> polinomio de fpx (objeto opaco) no nulo
# input: u --> polinomio irreducible mónico (objeto opaco) de grado >= 1
# output: multiplicidad de u como factor de f, es decir, el entero e >= 0
# mas grande tal que u^e | f
def multiplicidad_fpx(fpx, f, u):
    e = 0
    q, r = fpx.divmod(f, u)
    while fpx.es_cero(r):
        e += 1
        if fpx.es_cero(q):
            break
        f = q
        q, r = fpx.divmod(f, u)
    return e

# factorización de Cantor-Zassenhaus
# input: fpx --> anillo_fp_x (supondremos p impar)
# input: f --> polinomio de fpx (objeto opaco)
# output: [(f1,e1), ..., (fr,er)] donde f = c * f1^e1 * ... * fr^er es la
# factorización completa de f en irreducibles mónicos fi con multiplicidad
# ei >= 1 y los fi son distintos entre si y por supuesto c es el coeficiente
# principal de f
def fact_fpx(fpx: "cf.anillo_fp_x", f):
    g = sqfree_fact_fpx(fpx, f)
    h = didegr_fact_fpx(fpx, g)
    irreducibles = []
    for r in range(len(h)):
        if fpx.grado(h[r]) > 0:
            irreducibles += eqdegr_fact_fpx(fpx, r+1, h[r])
    factorizacion = []
    for u in irreducibles:
        e = multiplicidad_fpx(fpx, f, u)
        factorizacion += [(u,e)]
    return factorizacion

# esta linea es para añadir la función de factorización de Cantor-Zassenhaus
# como un método de la clase anillo_fp_x
cf.anillo_fp_x.factorizar = fact_fpx

# square-free factorization
# input: fqx --> anillo_fq_x
# input: f --> polinomio fabricado por fqx (objeto opaco) no nulo
# output: g = producto de los factores irreducibles mónicos de f, es decir,
# si f = c * f1^e1 * f2^e2 * ... * fr^er con los fi irreducibles mónicos
# distintos entre si, ei >= 1, c en fq, entonces g = f1 * f2 * ... * fr
def sqfree_fact_fqx(fqx, f):
    F = fqx._to_list(f)
    p = fqx.Q.p
    q_size = fqx.Q.fq

    # derivada (vectores)
    if len(F) <= 1:
        df = [fqx._zero_vec()]
    else:
        df = []
        for i in range(1, len(F)):
            k = i % p
            if k == 0:
                df.append(fqx._zero_vec())
            elif k == 1:
                df.append(F[i][:])
            else:
                df.append([ (k*c) % p for c in F[i] ])
    df = fqx._trim(df)

    def raiz_p_vec(vec):
        exp = q_size // p
        return fqx.Q._pow_vec(vec, exp)

    if fqx._gradoL(df) == -1:
        gL = []
        for i in range (0, len(F), p):
            gL.append(raiz_p_vec(F[i]))

        g_poly = fqx._from_list(fqx._trim(gL))
        return sqfree_fact_fqx(fqx, g_poly)

    cL = fqx._gcdL(F, df)
    qL, _ = fqx._divmodL(F, cL)
    gL = fqx._trim(qL)

    if fqx._gradoL(gL) > 0:
        lead_vec = gL[fqx._gradoL(gL)]
        if not fqx._is_uno_vec(lead_vec):
            inv_lead_vec = fqx.Q._pow_vec(lead_vec, q_size - 2)
            gL = [fqx.Q._mult_vec(c_vec, inv_lead_vec) for c_vec in gL]

    residuoL = [v[:] for v in F]
    while True:
        yL = fqx._gcdL(residuoL, gL)
        if fqx._gradoL(yL) <= 0:
            break
        residuoL, _ = fqx._divmodL(residuoL, yL)

    if fqx._gradoL(residuoL) > 0:
        h_coef = []
        for i in range(0, len(residuoL), p):
            h_coef.append(raiz_p_vec(residuoL[i]))

        h_poly = fqx._from_list(fqx._trim(h_coef))

        g2_poly = sqfree_fact_fqx(fqx, h_poly)

        g_poly = fqx._from_list(gL)

        # Unimos los factores de ambas partes evitando duplicados (LCM)
        producto = g_poly * g2_poly
        interseccion = fqx.gcd(g_poly, g2_poly)
        # producto / interseccion = LCM(g_poly, g2_poly)
        return fqx.div(producto, interseccion)
    
    else:
        
        return fqx._from_list(gL)



# distinct-degree factorization
# input: fqx --> anillo_fq_x
# input: g --> polinomio de fqx (objeto opaco) que es producto de factores
# irreducibles mónicos distintos cada uno con multiplicidad uno
# output: [h1, h2, ..., hr], donde hi = producto de los factores irreducibles
# mónicos de h de grado = i, el último hr debe ser no nulo y por supuesto
# g = h1 * h2 * ... * hr
def didegr_fact_fqx(fqx, g):
    q = fqx.Q.fq
    xL = [fqx._zero_vec(), fqx._one_vec()]
    gL = fqx._to_list(g)

    H_list = []
    H_current = fqx._modL(xL, gL)
    i = 1

    deg_g = fqx._gradoL(gL)

    while i<= deg_g and deg_g>0:
        H_current = fqx._pot_modL(H_current, q, gL)
        dL = fqx._gcdL(gL, fqx._sumaL(H_current, fqx._negL(xL)))

        if not fqx._is_unoL(dL):
            H_list.append(fqx._from_list(dL))
            qL, _ = fqx._divmodL(gL, dL)
            gL = qL
            deg_g = fqx._gradoL(gL)
            H_current = fqx._modL(H_current, gL) if deg_g >= 0 else [fqx._zero_vec()]

        else:
            H_list.append(fqx.uno())

        i += 1

        if i > deg_g // 2 and deg_g > 0:
            while len(H_list) < deg_g -1:
                H_list.append(fqx.uno())
            H_list.append(fqx._from_list(gL))
            break

    return H_list

# equal-degree factorization
# input: fqx --> anillo_fq_x (supondremos q impar)
# input: r --> int
# input: h --> polinomio de fqx (objeto opaco) que es producto de factores
# irreducibles mónicos distintos de grado r con multiplicidad uno
# output: [u1, ..., us], donde h = u1 * u2* ... * us y los ui son irreducibles
# mónicos de grado = r
def eqdegr_fact_fqx(fqx, r, h):
    q = fqx.Q.fq
    hL = fqx._to_list(h)
    deg_h = fqx._gradoL(hL)

    if deg_h <= r:
        if deg_h > 0:
            return [h]
        else:
            return []
    
    factors = []
    stack = [hL]

    while stack: 
        gL = stack.pop()
        deg_g = fqx._gradoL(gL)

        if deg_g == r:
            if deg_g > 0:
                factors.append(fqx._from_list(gL))
            continue

        if deg_g <= 0:
            continue
        
        found_split = False
        for _ in range(max(10, deg_g)):
            coefsv = [fqx._to_vec(fqx.Q.aleatorio()) for _ in range(deg_g)]
            coefsv.append(fqx._zero_vec())
            aL = fqx._trim(coefsv)
            aL = fqx._modL(aL, gL)

            if fqx._gradoL(aL) <=0:
                continue

            exp = (pow(q, r) -1)//2
            bL = fqx._pot_modL(aL, exp, gL)

            b_minus_1 = fqx._sumaL(bL, fqx._negL([fqx._one_vec()]))
            dL = fqx._gcdL(gL, b_minus_1)

            dr = fqx._gradoL(dL)

            if 0 < dr < deg_g:
                qL, _ = fqx._divmodL(gL, dL)
                stack.append(dL)
                stack.append(qL)
                found_split = True
                break

        if not found_split:
            if fqx._gradoL(gL) > 0:
                factors.append(fqx._from_list(gL))  

    return factors

# multiplicidad de factor irreducible mónico
# input: fqx --> anillo_fq_x
# input: f --> polinomio de fqx (objeto opaco) no nulo
# input: u --> polinomio irreducible mónico (objeto opaco) de grado >= 1
# output: multiplicidad de u como factor de f, es decir, el entero e >= 0
# mas grande tal que u^e | f
def multiplicidad_fqx(fqx, f, u):
    e = 0
    q, r = fqx.divmod(f, u)
    while fqx.es_cero(r):
        e += 1
        if fqx.es_cero(q):
            break
        f = q
        q, r = fqx.divmod(f, u)
    return e

# factorización de Cantor-Zassenhaus
# input: fqx --> anillo_fq_x (supondremos q impar)
# input: f --> polinomio de fqx (objeto opaco)
# output: [(f1,e1), ..., (fr,er)] donde f = c * f1^e1 * ... * fr^er es la
# factorización completa de f en irreducibles mónicos fi con multiplicidad
# ei >= 1 y los fi son distintos entre si y por supuesto c es el coeficiente
# principal de f
def fact_fqx(fqx: "cf.anillo_fp_x", f):                     # mantener esta implementación
    g = sqfree_fact_fqx(fqx, f)
    h = didegr_fact_fqx(fqx, g)
    irreducibles = []
    for r in range(len(h)):
        if fqx.grado(h[r]) > 0:
            irreducibles += eqdegr_fact_fqx(fqx, r+1, h[r])
    factorizacion = []
    for u in irreducibles:
        e = multiplicidad_fqx(fqx, f, u)
        factorizacion += [(u,e)]
    return factorizacion

# esta linea es para añadir la función de factorización de Cantor-Zassenhaus
# como un método de la clase anillo_fq_x
cf.anillo_fq_x.factorizar = fact_fqx
