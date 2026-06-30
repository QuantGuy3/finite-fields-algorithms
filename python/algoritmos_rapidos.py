# para esta actividad necesitamos la primera (cuerpos_finitos.py) completa
import unittest
import cuerpos_finitos as cf


"""
PRÁCTICA ALCP 3

Nombre y apellidos: Elia Torres Simón (eliatorr@ucm.es) 

Grupo: Ismael Amador, Ángela Ruiz, Miquel Sandonís, Elia Torres y Claudia Vicario
"""


# input: fpx -> anillo_fp_x
# input: f -> polinomio (objeto opaco creado por fpx)
# input: g -> polinomio (objeto opaco creado por fpx)
# output: f*g calculado usando el método de Karatsuba
def fp_x_mult_karatsuba(fpx, f, g):
    """
    Multiplica dos polinomios en Fp[x] usando Karatsuba.
    """
    # Si algún polinomio es cero, resultado es cero
    if fpx.es_cero(f) or fpx.es_cero(g):
        return fpx.cero()
    
    # Convertir a listas de coeficientes
    A = fpx._to_list(f)
    B = fpx._to_list(g)
    
    # Caso base: multiplicación directa para polinomios pequeños
    if len(A) == 1 or len(B) == 1:
        return fpx.mult(f, g)
    
    # Punto de división (mitad del tamaño mayor)
    m = max(len(A), len(B)) // 2
    
    # Dividir ambos polinomios en partes alta y baja
    def split_poly(P, split_point):
        if len(P) <= split_point:
            return P, [0]
        low = fpx._trim(P[:split_point])
        high = fpx._trim(P[split_point:])
        return low, high

    A_low, A_high = split_poly(A, m)  # A = A_low + A_high * x^m2
    B_low, B_high = split_poly(B, m)  # B = B_low + B_high * x^m2
    
    # Convertir partes a PolyK
    A_low_poly = fpx._from_list(A_low)
    A_high_poly = fpx._from_list(A_high)
    B_low_poly = fpx._from_list(B_low)
    B_high_poly = fpx._from_list(B_high)
    
    # --- ALGORITMO KARATSUBA (3 multiplicaciones en lugar de 4) ---
    
    # 1. Multiplicar partes bajas
    z0 = fp_x_mult_karatsuba(fpx, A_low_poly, B_low_poly)
    
    # 2. Multiplicar partes altas  
    z2 = fp_x_mult_karatsuba(fpx, A_high_poly, B_high_poly)
    
    # 3. Truco de Karatsuba: calcular término medio con 1 multiplicación
    # z1 = (A_low + A_high) * (B_low + B_high) - z0 - z2
    A_sum = fpx.suma(A_low_poly, A_high_poly)
    B_sum = fpx.suma(B_low_poly, B_high_poly)
    z1 = fp_x_mult_karatsuba(fpx, A_sum, B_sum)
    z1 = fpx.suma(z1, fpx.inv_adit(z0))  # - z0
    z1 = fpx.suma(z1, fpx.inv_adit(z2))  # - z2
    
    # Combinar resultados: z2*x^(2m2) + z1*x^m2 + z0
    z1_shifted = fpx._from_list([0] * m + fpx._to_list(z1))
    z2_shifted = fpx._from_list([0] * (2 * m) + fpx._to_list(z2))
    
    result = fpx.suma(z0, z1_shifted)
    result = fpx.suma(result, z2_shifted)
    
    return result

# añadimos esta función a la clase (sin sobreescribir la que ya teníamos)
cf.anillo_fp_x.mult_fast = fp_x_mult_karatsuba

# input: fp -> cuerpo_fp
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fp (primera columna de una
#    matriz de Toeplitz inferior T de nxn)
# input: b -> tupla de longitud n de elementos de fp (vector)
# output: T*b -> tupla de longitud n de elementos de fp (vector)
# se debe utilizar fp_x_mult_karatsuba internamente
def fp_toep_inf_vec(fp, n, a, b):
    fpx = cf.anillo_fp_x(fp)

    poly_a = fpx.elem_de_tuple(a)
    poly_b = fpx.elem_de_tuple(b)

    poly_c = fpx.mult_fast(poly_a, poly_b)

    all_coeffs = fpx.conv_a_tuple(poly_c)

    result = all_coeffs[:n]
    if len(result) < n:
        padding = (fp.cero(),) * (n-len(result))
        result += padding

    return result

# input: fp -> cuerpo_fp
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fp (primera fila de una
#    matriz de Toeplitz superior T de nxn)
# input: b -> tupla de longitud n de elementos de fp (vector)
# output: T*b -> tupla de longitud n de elementos de fp (vector)
# se debe utilizar fp_x_mult_karatsuba internamente
def fp_toep_sup_vec(fp, n, a, b):
    fpx = cf.anillo_fp_x(fp)
    # y = T_sup(a) * b
    # y_rev = T_inf(a) * b_rev
    b_rev = tuple(reversed(b))
    y_rev = fp_toep_inf_vec(fp, n, a, b_rev)
    return tuple(reversed(y_rev))


# input: fp -> cuerpo_fp
# input: n >= 1 (int)
# input: a -> tupla de longitud 2*n-1 de elementos de fp (primera fila de una
#    matriz de Toeplitz completa T de nxn seguida de la primera columna
#    excepto el elemento de la esquina)
# input: b -> tupla de longitud n de elementos de fp (vector)
# output: T*b -> tupla de longitud n de elementos de fp (vector)
# se debe utilizar fp_x_mult_karatsuba internamente
def fp_toep_vec(fp, n, a, b):
    fpx = cf.anillo_fp_x(fp)

    a_neg = tuple(reversed(a[1:n]))
    a_0 = (a[0],)
    a_pos = a[n:]
    poly_a_coeffs = a_neg + a_0 + a_pos

    poly_a = fpx.elem_de_tuple(poly_a_coeffs)
    poly_b = fpx.elem_de_tuple(b)

    poly_c = fpx.mult_fast(poly_a, poly_b)

    all_coeffs = fpx.conv_a_tuple(poly_c)

    len_necesaria = 2 * n - 1
    if len(all_coeffs) < len_necesaria:
        padding = (fp.cero(),) * (len_necesaria - len(all_coeffs))
        all_coeffs += padding

    result_coeffs = all_coeffs[n-1:2*n-1]
    return result_coeffs

# input: fp -> cuerpo_fp
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fp (primera columna de una
#    matriz de Toeplitz inferior T de nxn)... suponemos a[0] != 0
# output: primera columna de T^(-1) -> tupla de longitud n de elementos de
#    fp (vector)
# utilizar un método recursivo que "divida el problema a la mitad"
# recordar que T^(-1) es también una matriz de Toeplitz inferior
def fp_toep_inf_inv(fp, n, a):
    if n == 1:
        if fp.es_cero(a[0]):
            raise ZeroDivisionError("El elemento a[0] no puede ser cero.")
        return (fp.inv_mult(a[0]),)
    
    k = (n+1)//2  

    a_k = a[:k]
    x_k_coeffs = fp_toep_inf_inv(fp, k, a_k)

    fpx = cf.anillo_fp_x(fp)

    poly_Xk = fpx.elem_de_tuple(x_k_coeffs)
    poly_An = fpx.elem_de_tuple(a)

    poly_T = fpx.mult_fast(poly_An, poly_Xk)

    poly_2 = fpx.elem_de_int(2)
    poly_T_neg = fpx.inv_adit(poly_T)
    poly_R = fpx.suma(poly_T_neg, poly_2)

    poly_Xn_full = fpx.mult_fast(poly_R, poly_Xk)

    all_coeffs = fpx.conv_a_tuple(poly_Xn_full)

    result = all_coeffs[:n]
    if len(result) < n:
        padding = (fp.cero(),) * (n-len(result))
        result += padding

    return result

# input: fp -> cuerpo_fp
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fp (primera fila de una
#    matriz de Toeplitz superior T de nxn)... suponemos a[0] != 0
# output: primera fila de T^(-1) -> tupla de longitud n de elementos de
#    fp (vector)
# utilizar un método recursivo que "divida el problema a la mitad"
# recordar que T^(-1) es también una matriz de Toeplitz superior
def fp_toep_sup_inv(fp, n, a):
    return fp_toep_inf_inv(fp, n, a)

# input: fpx -> anillo_fp_x
# input: f -> polinomio (objeto opaco creado por fpx)
# input: g -> polinomio no nulo (objeto opaco creado por fpx)
# output: q -> cociente
# output: r -> resto
# se cumple que f = g*q+r, r=0 o deg(r)<deg(g)
# reformular el problema en términos de matrices de Toeplitz y luego usar
# las funciones de arriba para obtener q y r
def fp_x_divmod(fpx, f, g):
    # Coeficientes en Z_p
    F = fpx._to_list(f)
    G = fpx._to_list(g)
    m = fpx._gradoL(F)
    n = fpx._gradoL(G)

    if n == -1:
        raise ZeroDivisionError("División por 0 en Fp[x]")

    # deg(f) < deg(g) ⇒ q = 0, r = f
    if m < n:
        return fpx.cero(), f

    p = fpx.p
    d = m - n   # grado de q

    # f_rev[i] = coef de x^{m-i} en f
    f_rev = [0] * (m + 1)
    for i, ci in enumerate(F):
        f_rev[m - i] = ci

    # g_rev[j] = coef de x^{n-j} en g
    g_rev = [0] * (n + 1)
    for j, cj in enumerate(G):
        g_rev[n - j] = cj

    # f_hat = primeros d+1 coeficientes de f_rev
    f_hat = f_rev[:d + 1]
    q_hat = [0] * (d + 1)

    g0 = g_rev[0]             # coeficiente líder de g
    inv_g0 = fpx.K._inv(g0)   # inverso en Fp

    # Sistema triangular: para i = 0..d
    # g0 q_hat[i] + sum_{j=1..min(i,n)} g_rev[j] q_hat[i-j] = f_hat[i]
    for i in range(d + 1):
        s = 0
        j_max = min(i, n)
        for j in range(1, j_max + 1):
            s = (s + g_rev[j] * q_hat[i - j]) % p
        q_hat[i] = ((f_hat[i] - s) * inv_g0) % p

    # Recuperar q(x): q_k = q_hat[d-k]
    q_coeffs = [0] * (d + 1)
    for k in range(d + 1):
        q_coeffs[k] = q_hat[d - k]

    # r(x) = f(x) - g(x)*q(x)
    gq = [0] * (m + 1)      # deg(gq) ≤ m
    for i, ai in enumerate(q_coeffs):
        for j, bj in enumerate(G):
            gq[i + j] = (gq[i + j] + ai * bj) % p

    # ajustamos longitudes (deberían coincidir ya)
    if len(F) < len(gq):
        F += [0] * (len(gq) - len(F))
    if len(gq) < len(F):
        gq += [0] * (len(F) - len(gq))

    r_coeffs = [(fi - gi) % p for fi, gi in zip(F, gq)]

    # normalizar
    q_coeffs = fpx._trim(q_coeffs)
    r_coeffs = fpx._trim(r_coeffs)

    q_poly = fpx._from_list(q_coeffs)
    r_poly = fpx._from_list(r_coeffs)
    return q_poly, r_poly



# añadimos esta función a la clase (sin sobreescribir la que ya teníamos)
cf.anillo_fp_x.divmod_fast = fp_x_divmod

# input: fp -> cuerpo_fp
# input: g -> elemento del grupo multiplicativo fp* de orden n (objeto opaco)
# input: k >= 0 tal que n = 2**k divide a p-1
# input: a -> tupla de longitud n de elementos de fp
# output: DFT_{n,g}(a) -> tupla de longitud n de elementos de fp
# utilizar el algoritmo de Cooley-Tuckey
def fp_fft(fp, g, k, a):
    n = 2**k
    if len(a) != n:
        raise ValueError("fp_fft: len(a) debe ser 2**k")
    if n == 1: 
        return tuple([a[0]])

    a_pares = a[0::2]
    a_impares  = a[1::2]
    
    g2 = fp.pot(g, 2)  # FFT recursiva con raíz g^2
    P = fp_fft(fp, g2, k - 1, a_pares)
    I = fp_fft(fp, g2, k - 1, a_impares)
    
    A = [fp.cero() for _ in range(n)]
    w = fp.uno()  
    for j in range(n // 2):
         t = fp.mult(w, I[j])
         A[j] = fp.suma(P[j], t)
         A[j + n // 2] = fp.suma(P[j], fp.inv_adit(t))
         w = fp.mult(w, g)
    return tuple(A)

    
# input: fp -> cuerpo_fp
# input: g -> elemento del grupo multiplicativo fp* de orden n (objeto opaco)
# input: k >= 0 tal que n = 2**k divide a p-1
# input: a -> tupla de longitud n de elementos de fp
# output: IDFT_{n,g}(a) -> tupla de longitud n de elementos de fp
# recordar que IDFT_{n,g} = n^(-1) * DFT_{n,g^(-1)}
def fp_ifft(fp, g, k, a):
    n = 2 ** k
    if len(a) != n:
        raise ValueError("fp_ifft: len(a) debe ser 2**k")
        
    g_inv = fp.inv_mult(g)
    A = fp_fft(fp, g_inv, k, a) # FFT con la raíz inversa
    n_inv = fp.inv_mult(fp.elem(n))
    return tuple(fp.mult(n_inv, x) for x in A)

# input: fqx -> anillo_fq_x
# input: f -> polinomio (objeto opaco creado por fqx)
# input: g -> polinomio (objeto opaco creado por fqx)
# output: f*g calculado usando el método de Karatsuba
def fq_x_mult_karatsuba(fqx, f, g):
    """
    Multiplica dos polinomios en Fq[x] usando Karatsuba.
    """
    # Si algún polinomio es cero, resultado es cero
    if fqx.es_cero(f) or fqx.es_cero(g):
        return fqx.cero()
    
    # Convertir a listas de vectores (cada vector representa un elemento de Fq)
    A = fqx._to_list(f)
    B = fqx._to_list(g)
    
    # Caso base: multiplicación directa para polinomios pequeños
    if len(A) == 1 or len(B) == 1:
        return fqx.mult(f, g)
    
    # Punto de división (mitad del tamaño mayor)
    m = max(len(A), len(B)) // 2
    
    # Dividir ambos polinomios en partes alta y baja
    def split_poly(P, split_point):
        if len(P) <= split_point:
            return fqx._trim(P), [fqx._zero_vec()]
        low  = fqx._trim(P[:split_point])
        high = fqx._trim(P[split_point:])
        return low, high
    
    A_low, A_high = split_poly(A, m)  # A = A_low + A_high * x^m2
    B_low, B_high = split_poly(B, m)  # B = B_low + B_high * x^m2
    
    # Convertir partes a PolyK
    A_low_poly = fqx._from_list(A_low)
    A_high_poly = fqx._from_list(A_high)
    B_low_poly = fqx._from_list(B_low)
    B_high_poly = fqx._from_list(B_high)
    
    # --- ALGORITMO KARATSUBA (3 multiplicaciones en lugar de 4) ---
    
    # 1. Multiplicar partes bajas
    z0 = fq_x_mult_karatsuba(fqx, A_low_poly, B_low_poly)
    
    # 2. Multiplicar partes altas  
    z2 = fq_x_mult_karatsuba(fqx, A_high_poly, B_high_poly)
    
    # 3. Truco de Karatsuba: calcular término medio con 1 multiplicación
    # z1 = (A_low + A_high) * (B_low + B_high) - z0 - z2
    A_sum = fqx.suma(A_low_poly, A_high_poly)
    B_sum = fqx.suma(B_low_poly, B_high_poly)
    z1 = fq_x_mult_karatsuba(fqx, A_sum, B_sum)
    z1 = fqx.suma(z1, fqx.inv_adit(z0))  # - z0
    z1 = fqx.suma(z1, fqx.inv_adit(z2))  # - z2
    
    # Combinar resultados: z2*x^(2m2) + z1*x^m2 + z0
    z1_shifted = fqx._from_list([fqx._zero_vec()] * m + fqx._to_list(z1))
    z2_shifted = fqx._from_list([fqx._zero_vec()] * (2 * m) + fqx._to_list(z2))
    
    result = fqx.suma(z0, z1_shifted)
    result = fqx.suma(result, z2_shifted)
    
    return result

# añadimos esta función a la clase (sin sobreescribir la que ya teníamos)
cf.anillo_fq_x.mult_fast = fq_x_mult_karatsuba

# input: fq -> cuerpo_fq
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fq (primera columna de una
#    matriz de Toeplitz inferior T de nxn)
# input: b -> tupla de longitud n de elementos de fq (vector)
# output: T*b -> tupla de longitud n de elementos de fq (vector)
# se debe utilizar fq_x_mult_karatsuba internamente
def fq_toep_inf_vec(fq, n, a, b):
    fqx = cf.anillo_fq_x(fq)

    poly_a = fqx.elem_de_tuple(a)
    poly_b = fqx.elem_de_tuple(b)

    poly_c = fqx.mult_fast(poly_a, poly_b)

    all_coeffs = fqx.conv_a_tuple(poly_c)

    result = all_coeffs[:n]
    if len(result) < n:
        padding = (fq.cero(),) * (n-len(result))
        result += padding

    return result

# input: fq -> cuerpo_fq
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fq (primera fila de una
#    matriz de Toeplitz superior T de nxn)
# input: b -> tupla de longitud n de elementos de fq (vector)
# output: T*b -> tupla de longitud n de elementos de fq (vector)
# se debe utilizar fq_x_mult_karatsuba internamente
def fq_toep_sup_vec(fq, n, a, b):
    fqx = cf.anillo_fq_x(fq)
    b_rev = tuple(reversed(b))
    y_rev = fq_toep_inf_vec(fq, n, a, b_rev)
    return tuple(reversed(y_rev))


# input: fq -> cuerpo_fq
# input: n >= 1 (int)
# input: a -> tupla de longitud 2*n-1 de elementos de fq (primera fila de una
#    matriz de Toeplitz completa T de nxn seguida de la primera columna
#    excepto el elemento de la esquina)
# input: b -> tupla de longitud n de elementos de fq (vector)
# output: T*b -> tupla de longitud n de elementos de fq (vector)
# se debe utilizar fq_x_mult_karatsuba internamente
def fq_toep_vec(fq, n, a, b):
    fqx = cf.anillo_fq_x(fq)

    a_neg = tuple(reversed(a[1:n]))
    a_0 = (a[0],)
    a_pos = a[n:]
    poly_a_coeffs = a_neg + a_0 + a_pos

    poly_a = fqx.elem_de_tuple(poly_a_coeffs)
    poly_b = fqx.elem_de_tuple(b)

    poly_c = fqx.mult_fast(poly_a, poly_b)

    all_coeffs = fqx.conv_a_tuple(poly_c)

    len_necesaria = 2 * n - 1
    if len(all_coeffs) < len_necesaria:
        padding = (fq.cero(),) * (len_necesaria - len(all_coeffs))
        all_coeffs += padding

    result = all_coeffs[n-1 : 2*n-1]
    return result

    # FALTA

# input: fq -> cuerpo_fq
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fq (primera columna de una
#    matriz de Toeplitz inferior T de nxn)... suponemos a[0] != 0
# output: primera columna de T^(-1) -> tupla de longitud n de elementos de
#    fq (vector)
# utilizar un método recursivo que "divida el problema a la mitad"
# recordar que T^(-1) es también una matriz de Toeplitz inferior
def fq_toep_inf_inv(fq, n, a):
    if n == 1:
        if fq.es_cero(a[0]):
            raise ZeroDivisionError("Matriz singular, a[0] es cero")
        # columna de tamaño 1
        return (fq.inv_mult(a[0]),)

    k = (n + 1) // 2

    a_k = a[:k]
    x_k_coeffs = fq_toep_inf_inv(fq, k, a_k)

    fqx = cf.anillo_fq_x(fq)

    poly_Xk = fqx.elem_de_tuple(x_k_coeffs)
    poly_An = fqx.elem_de_tuple(a)

    # T(Xk) ~ A_n * X_k en la convolución
    poly_T = fqx.mult_fast(poly_An, poly_Xk)

    poly_T_neg = fqx.inv_adit(poly_T)

    # 2 en Fq: uno + uno (en característica 2 esto da 0, como debe ser)
    two = fq.suma(fq.uno(), fq.uno())
    poly_2 = fqx.elem_de_tuple((two,))

    # R = 2 - T(Xk)
    poly_R = fqx.suma(poly_T_neg, poly_2)

    # X_n = R * X_k
    poly_Xn_full = fqx.mult_fast(poly_R, poly_Xk)

    all_coeffs = fqx.conv_a_tuple(poly_Xn_full)

    result_coeffs = all_coeffs[:n]
    if len(result_coeffs) < n:
        padding = (fq.cero(),) * (n - len(result_coeffs))
        result_coeffs += padding

    return result_coeffs



# input: fq -> cuerpo_fq
# input: n >= 1 (int)
# input: a -> tupla de longitud n de elementos de fq (primera fila de una
#    matriz de Toeplitz superior T de nxn)... suponemos a[0] != 0
# output: primera fila de T^(-1) -> tupla de longitud n de elementos de
#    fq (vector)
# utilizar un método recursivo que "divida el problema a la mitad"
# recordar que T^(-1) es también una matriz de Toeplitz superior
def fq_toep_sup_inv(fq, n, a):
    return fq_toep_inf_inv(fq, n, a)

# input: fqx -> anillo_fq_x
# input: f -> polinomio (objeto opaco creado por fqx)
# input: g -> polinomio no nulo (objeto opaco creado por fqx)
# output: q -> cociente
# output: r -> resto
# se cumple que f = g*q+r, r=0 o deg(r)<deg(g)
# reformular el problema en términos de matrices de Toeplitz y luego usar
# las funciones de arriba para obtener q y r
def fq_x_divmod(fqx, f, g):
    fq = fqx.Q  # cuerpo Fq

    # Coeficientes como ElementoFq
    F = list(fqx.conv_a_tuple(f))
    G = list(fqx.conv_a_tuple(g))

    # quitar ceros de cola
    while len(F) > 1 and fq.es_cero(F[-1]):
        F.pop()
    while len(G) > 1 and fq.es_cero(G[-1]):
        G.pop()

    m = len(F) - 1
    n = len(G) - 1

    if n == -1:
        raise ZeroDivisionError("División por 0 en Fq[x]")

    if m < n:
        return fqx.cero(), f   # q = 0, r = f

    d = m - n  # grado de q

    # f_rev[i] = coef de x^{m-i} en f
    f_rev = [fq.cero()] * (m + 1)
    for i, ci in enumerate(F):
        f_rev[m - i] = ci

    # g_rev[j] = coef de x^{n-j} en g
    g_rev = [fq.cero()] * (n + 1)
    for j, cj in enumerate(G):
        g_rev[n - j] = cj

    f_hat = f_rev[:d + 1]
    q_hat = [fq.cero() for _ in range(d + 1)]

    g0 = g_rev[0]
    g0_inv = fq.inv_mult(g0)

    # g0 q_hat[i] + sum_{j=1..min(i,n)} g_rev[j] q_hat[i-j] = f_hat[i]
    for i in range(d + 1):
        s = fq.cero()
        j_max = min(i, n)
        for j in range(1, j_max + 1):
            s = fq.suma(s, fq.mult(g_rev[j], q_hat[i - j]))
        q_hat[i] = fq.mult(fq.suma(f_hat[i], fq.inv_adit(s)), g0_inv)

    # q(x) a partir de q_hat
    q_coeffs = [fq.cero() for _ in range(d + 1)]
    for k in range(d + 1):
        q_coeffs[k] = q_hat[d - k]

    q_poly = fqx.elem_de_tuple(tuple(q_coeffs))

    # resto r = f - g*q
    gq = fqx.mult(g, q_poly)           # g*q
    r_poly = fqx.suma(f, fqx.inv_adit(gq))

    # normalizar coeficientes del resto
    r_coeffs = list(fqx.conv_a_tuple(r_poly))
    while len(r_coeffs) > 1 and fq.es_cero(r_coeffs[-1]):
        r_coeffs.pop()
    r_poly = fqx.elem_de_tuple(tuple(r_coeffs))

    return q_poly, r_poly



# añadimos esta función a la clase (sin sobreescribir la que ya teníamos)
cf.anillo_fq_x.divmod_fast = fq_x_divmod

# input: fq -> cuerpo_fq
# input: g -> elemento del grupo multiplicativo fq* de orden n (objeto opaco)
# input: k >= 0 tal que n = 2**k divide a q-1
# input: a -> tupla de longitud n de elementos de fq
# output: DFT_{n,g}(a) -> tupla de longitud n de elementos de fq
# utilizar el algoritmo de Cooley-Tuckey
def fq_fft(fq, g, k, a):
    n = 2**k
    if len(a) != n:
        raise ValueError("fq_fft: len(a) debe ser 2**k")
    if n == 1: 
        return tuple([a[0]])

    a_pares = a[0::2]
    a_impares  = a[1::2]
    
    g2 = fq.pot(g, 2)
    P = fq_fft(fq, g2, k - 1, a_pares)
    I = fq_fft(fq, g2, k - 1, a_impares)
    
    A = [fq.cero() for _ in range(n)]
    w = fq.uno()  # g^0
    for j in range(n // 2):
        t = fq.mult(w, I[j])
        A[j] = fq.suma(P[j], t)
        A[j + n // 2] = fq.suma(P[j], fq.inv_adit(t))
        w = fq.mult(w, g)
    
    return tuple(A)

# input: fq -> cuerpo_fq
# input: g -> elemento del grupo multiplicativo fq* de orden n (objeto opaco)
# input: k >= 0 tal que n = 2**k divide a p-1
# input: a -> tupla de longitud n de elementos de fq
# output: IDFT_{n,g}(a) -> tupla de longitud n de elementos de fq
# recordar que IDFT_{n,g} = n^(-1) * DFT_{n,g^(-1)}
def fq_ifft(fq, g, k, a):
    n = 2**k
    if len(a) != n:
        raise ValueError("fq_ifft: len(a) debe ser 2**k")
        
    g_inv = fq.inv_mult(g)
    A_fft = fq_fft(fq, g_inv, k, a) # FFT con la raíz inversa
    n_inv = fq.inv_mult(fq.elem_de_int(n))
    A_final = [fq.mult(n_inv, x) for x in A_fft]
    
    return tuple(A_final)

























