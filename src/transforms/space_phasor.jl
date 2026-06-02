# =============================================================================
# Fasor espacial de corriente (fmm) para máquinas de campo giratorio
#
# Derivación (§1.6):
#   El devanado trifásico simétrico produce distribuciones sinusoidales de fmm
#   en el entrehierro. Cada fase se representa como un fasor en el plano complejo
#   cuyo eje real coincide con el eje magnético de la fase a.
#
#   Suma vectorial de los tres fasores de fase:
#       is_raw = ia + a·ib + a²·ic       con  a = e^(j2π/3)
#
#   El factor k = 2/3 se elige para que la proyección sobre el eje de la fase a
#   reproduzca el valor instantáneo de la corriente en esa fase:
#       Re{is} = (2/3)·Re{ia + a·ib + a²·ic} = ia      (cuando ia+ib+ic = 0)
#
#   Definición final (ec. 1.6.4):
#       is = (2/3)·(ia + a·ib + a²·ic)
#
# Transformación de coordenadas (ec. 1.6.7):
#   El fasor visto desde un sistema girante a velocidad dγ/dt se obtiene
#   multiplicando por e^(-jγ):
#       is,r = is · e^(-jγ)
#
# Componentes reales (ec. 1.6.9 / 1.6.13):
#   Proyectando sobre los ejes del sistema girante (d,q):
#       id + j·iq = is,r = is · e^(-jγ)
#
# Componente de secuencia cero (ec. 1.6.12):
#   No contribuye al campo fundamental; debe especificarse por separado:
#       i0 = (1/3)·(ia + ib + ic)
# =============================================================================

const _a  = exp(im * 2π/3)   # operador de desfase 120°
const _a2 = exp(im * 4π/3)   # operador de desfase 240°

# -----------------------------------------------------------------------------
# Construcción del fasor espacial desde magnitudes de fase (abc)
# Ec. (1.6.4): is = (2/3)·(ia + a·ib + a²·ic)
# -----------------------------------------------------------------------------

"""
    spatial_phasor(ia, ib, ic) -> ComplexF64

Fasor espacial de corriente en el sistema de referencia fijo al estator
(eje real alineado con el eje magnético de la fase a).

La proyección sobre el eje real reproduce el valor instantáneo de ia:
    Re{is} = ia    (válido cuando ia + ib + ic = 0)
"""
function spatial_phasor(ia::Float64, ib::Float64, ic::Float64)
    return (2/3) * (ia + _a * ib + _a2 * ic)
end

# -----------------------------------------------------------------------------
# Transformación al sistema de referencia girante (eje d–q)
# Ec. (1.6.7): is,r = is · e^(-jγ)
# γ es el ángulo eléctrico entre el eje d del sistema girante y el eje de la fase a.
# -----------------------------------------------------------------------------

"""
    phasor_to_rotating(is, γ) -> ComplexF64

Transforma el fasor espacial `is` (en coordenadas fijas al estator)
al sistema de referencia girante cuyo eje d forma un ángulo γ con el eje α.

    is,r = is · e^(-jγ)
    id   = Re{is,r}
    iq   = Im{is,r}

Para el sistema fijo al rotor: γ = ángulo eléctrico del rotor.
Para el sistema sincrónico:    γ = ωs·t + γ₀.
"""
function phasor_to_rotating(is::ComplexF64, γ::Float64)
    return is * exp(-im * γ)
end

"""
    phasor_from_rotating(is_r, γ) -> ComplexF64

Transformación inversa: desde coordenadas girantes (d,q) a fijas al estator.

    is = is,r · e^(jγ)
"""
function phasor_from_rotating(is_r::ComplexF64, γ::Float64)
    return is_r * exp(im * γ)
end

# -----------------------------------------------------------------------------
# Componente de secuencia cero
# Ec. (1.6.12): i0 = (1/3)·(ia + ib + ic)
# No contribuye al campo fundamental en el entrehierro.
# Para devanado en estrella sin neutro: i0 = 0.
# -----------------------------------------------------------------------------

"""
    zero_sequence(ia, ib, ic) -> Float64

Componente de secuencia cero. Es nula para conexión en estrella sin neutro.
"""
zero_sequence(ia::Float64, ib::Float64, ic::Float64) = (ia + ib + ic) / 3

# -----------------------------------------------------------------------------
# Reconstrucción de las magnitudes de fase desde el fasor espacial
# Ec. (1.7.4)–(1.7.6): va = Re{vs},  vb = Re{a²·vs},  vc = Re{a·vs}
# Se requiere adicionalmente la componente de secuencia cero i0.
# -----------------------------------------------------------------------------

"""
    phasor_to_abc(is, i0) -> (ia, ib, ic)

Reconstruye las corrientes de fase a partir del fasor espacial y la componente
de secuencia cero (necesaria para que la transformación sea biunívoca).
"""
function phasor_to_abc(is::ComplexF64, i0::Float64)
    ia = real(is)           + i0
    ib = real(_a2 * is)     + i0
    ic = real(_a  * is)     + i0
    return ia, ib, ic
end

# -----------------------------------------------------------------------------
# Transformación directa abc → αβ (Clarke) expresada en términos del fasor
# El fasor espacial en coordenadas estacionarias equivale a:
#       isα = Re{is} = (2/3)·(ia - ib/2 - ic/2)
#       isβ = Im{is} = (2/3)·(√3/2·ib - √3/2·ic)
# que es exactamente la transformada de Clarke (factor 2/3, no normalizada).
# -----------------------------------------------------------------------------

"""
    phasor_to_alphabeta(is) -> (iα, iβ)

Extrae las componentes αβ del fasor espacial en coordenadas estacionarias.
Equivale a la transformada de Clarke con factor 2/3.
"""
phasor_to_alphabeta(is::ComplexF64) = (real(is), imag(is))

"""
    alphabeta_to_phasor(iα, iβ) -> ComplexF64

Construye el fasor espacial desde las componentes αβ.
"""
alphabeta_to_phasor(iα::Float64, iβ::Float64) = complex(iα, iβ)

# -----------------------------------------------------------------------------
# Transformación completa abc → dq (Park) vía fasor espacial
# Combina spatial_phasor + phasor_to_rotating en un solo paso.
# -----------------------------------------------------------------------------

"""
    abc_to_dq(ia, ib, ic, γ) -> (id, iq, i0)

Transformación de Park completa: abc → dq0 usando el fasor espacial.
    is   = (2/3)·(ia + a·ib + a²·ic)
    is,r = is · e^(-jγ)
    id   = Re{is,r},  iq = Im{is,r},  i0 = (ia+ib+ic)/3
"""
function abc_to_dq(ia::Float64, ib::Float64, ic::Float64, γ::Float64)
    is   = spatial_phasor(ia, ib, ic)
    is_r = phasor_to_rotating(is, γ)
    i0   = zero_sequence(ia, ib, ic)
    return real(is_r), imag(is_r), i0
end

"""
    dq_to_abc(id, iq, i0, γ) -> (ia, ib, ic)

Transformación de Park inversa: dq0 → abc usando el fasor espacial.
"""
function dq_to_abc(id::Float64, iq::Float64, i0::Float64, γ::Float64)
    is_r = complex(id, iq)
    is   = phasor_from_rotating(is_r, γ)
    return phasor_to_abc(is, i0)
end
