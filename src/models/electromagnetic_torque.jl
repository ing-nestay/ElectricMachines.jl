# =============================================================================
# Momento electromagnético — Fuerzas de Lorentz (§1.8)
#
# Derivación:
#   Se reemplaza el devanado de armadura por una capa de corriente de densidad
#   lineal a(x) [A/m] sobre la superficie interior del estator.
#   La inducción resultante en el entrehierro sea b(x).
#
#   Fuerza tangencial sobre un elemento diferencial (ec. 1.8.1–1.8.2):
#       df  = −b(x)·a(x)·l·R·dx
#       dTs = −R²·l·b(x)·a(x)·dx
#
#   Momento sobre el estator (ec. 1.8.3):
#       Ts = −R²·l · ∫₀²π b(x)·a(x) dx
#
#   Momento sobre el rotor (tercera ley de Newton, ec. 1.8.4):
#       T  = −Ts = R²·l · ∫₀²π b(x)·a(x) dx
#
#   Evaluando la integral para la fundamental de la onda de inducción
#   y expresando el resultado en términos de las variables dq (ec. 1.8.x):
#
#       T = (3/2)·p·(ψ1d·i1q − ψ1q·i1d)
#
#   Forma equivalente en términos del fasor espacial (ec. 1.42 del modelo MAS):
#       T = −(3/2)·p·Im{ψs · conj(is)}
#
#   Ecuación de movimiento (ec. 1.8, §1.8.1):
#       J·dω/dt = p·(T − Tm)
#       dγ/dt   = ω
# =============================================================================

# -----------------------------------------------------------------------------
# Torque en coordenadas dq
# Resultado de evaluar la integral de Lorentz con la fundamental de la fmm.
# El factor 3/2 proviene del promedio sobre las tres fases simétricas.
# -----------------------------------------------------------------------------

"""
    torque_dq(ψ1d, ψ1q, i1d, i1q, p) -> Float64

Momento electromagnético en coordenadas dq (ec. 1.8):
    T = (3/2)·p·(ψ1d·i1q − ψ1q·i1d)

Válido para máquina sincrónica y asincrónica en el marco dq.
El signo es positivo para funcionamiento motor (T > 0 → aceleración).
"""
function torque_dq(ψ1d::Float64, ψ1q::Float64,
                   i1d::Float64, i1q::Float64,
                   p::Int)
    return (3/2) * p * (ψ1d * i1q - ψ1q * i1d)
end

# -----------------------------------------------------------------------------
# Torque en términos del fasor espacial
# Equivalente a la expresión anterior pero usando variables complejas.
# Útil para el modelo de la MAS en coordenadas sincrónicas.
# -----------------------------------------------------------------------------

"""
    torque_phasor(ψs, is, p) -> Float64

Momento electromagnético en términos del fasor espacial (ec. 1.42):
    T = −(3/2)·p·Im{ψs · conj(is)}

Equivalente a torque_dq cuando ψs = ψ1d + j·ψ1q e is = i1d + j·i1q.
"""
function torque_phasor(ψs::ComplexF64, is::ComplexF64, p::Int)
    return -(3/2) * p * imag(ψs * conj(is))
end

# -----------------------------------------------------------------------------
# Ecuación de movimiento mecánica
# Ec. (1.8, §1.8.1): J·dω/dt = p·(T − Tm)
# Tm > 0 → carga resistente (motor); Tm < 0 → par motor externo (generador).
# -----------------------------------------------------------------------------

"""
    mechanical_ode(T, Tm, ω, J, p) -> (dω, dγ)

Derivadas del estado mecánico:
    dω/dt = (p/J)·(T − Tm)      [rad_elec/s²]
    dγ/dt = ω                   [rad_elec/s]

Entradas:
- T  : torque electromagnético [N·m]
- Tm : torque mecánico externo  [N·m]  (positivo = carga resistente)
- ω  : velocidad eléctrica      [rad_elec/s]
- J  : inercia total            [kg·m²]
- p  : pares de polos
"""
function mechanical_ode(T::Float64, Tm::Float64,
                         ω::Float64, J::Float64, p::Int)
    dω = (p / J) * (T - Tm)
    dγ = ω
    return dω, dγ
end

# -----------------------------------------------------------------------------
# Potencia electromagnética
# Se obtiene del producto torque × velocidad mecánica.
# ωm = ω/p  →  Pem = T·ωm = T·ω/p
# -----------------------------------------------------------------------------

"""
    electromagnetic_power(T, ω, p) -> Float64

Potencia electromagnética [W]:
    Pem = T · ω / p
donde ω es la velocidad angular eléctrica y p el número de pares de polos.
"""
electromagnetic_power(T::Float64, ω::Float64, p::Int) = T * ω / p
