# =============================================================================
# Modelos Motor/Generador DC
# Versión: Excitación Independiente (Separada)
# =============================================================================
# Estados:
#   x[1] = if  → corriente de campo [A]
#   x[2] = ia  → corriente de armadura [A]
#   x[3] = ω   → velocidad angular [rad/s]
#
# Entradas:
#   Vf  → tensión de campo [V]
#   Va  → tensión de armadura [V]
#   TL  → torque de carga [N·m]  (positivo = carga motora)
# =============================================================================

struct DCMotorIndParams
    Rf::Float64    # Resistencia de campo [Ω]
    Lf::Float64    # Inductancia de campo [H]
    Ra::Float64    # Resistencia de armadura [Ω]
    La::Float64    # Inductancia de armadura [H]
    Ke::Float64    # Constante de máquina [V·s/A·rad] (flujo por corriente de campo)
    J::Float64     # Inercia total del rotor [kg·m²]
    B::Float64     # Fricción viscosa [N·m·s/rad]
end

"""
    dc_ind_ode!(dx, x, p, t, Vf, Va, TL)

Ecuaciones de estado del motor DC de excitación independiente.

Estados x = [if, ia, ω]:
- x[1] = if  corriente de campo [A]
- x[2] = ia  corriente de armadura [A]
- x[3] = ω   velocidad angular [rad/s]
"""
function dc_ind_ode!(dx, x, p::DCMotorIndParams, t, Vf, Va, TL)
    if_curr, ia, ω = x[1], x[2], x[3]

    # Fuerza contra-electromotriz (back-EMF)
    E = p.Ke * if_curr * ω

    # Torque electromagnético
    Te = p.Ke * if_curr * ia

    # Ecuaciones diferenciales
    dx[1] = (Vf - p.Rf * if_curr) / p.Lf          # circuito de campo
    dx[2] = (Va - p.Ra * ia - E) / p.La            # circuito de armadura
    dx[3] = (Te - p.B * ω - TL) / p.J              # mecánica
end

"""
    dc_ind_steady_state(p, Vf, Va, TL) -> (if_ss, ia_ss, ω_ss, Te_ss)

Punto de operación en estado estacionario (todas las derivadas = 0).
"""
function dc_ind_steady_state(p::DCMotorIndParams, Vf::Float64, Va::Float64, TL::Float64)
    # Campo: Lf·dif/dt = 0 → if_ss = Vf/Rf
    if_ss = Vf / p.Rf

    # Sistema armadura + mecánica en SS:
    # Ra·ia = Va - Ke·if·ω
    # Ke·if·ia = B·ω + TL
    # Resolviendo: ω_ss = (Ke·if·Va - Ra·TL) / (Ke²·if² + Ra·B)
    KΦ = p.Ke * if_ss
    ω_ss = (KΦ * Va - p.Ra * TL) / (KΦ^2 + p.Ra * p.B)
    ia_ss = (Va - KΦ * ω_ss) / p.Ra
    Te_ss = KΦ * ia_ss

    return if_ss, ia_ss, ω_ss, Te_ss
end

"""
    dc_ind_torque_speed(p, Vf, Va; TL_range) -> (ω_vec, Te_vec)

Genera la curva par-velocidad para excitación independiente.
TL_range: rango de torques de carga a evaluar.
"""
function dc_ind_torque_speed(p::DCMotorIndParams, Vf::Float64, Va::Float64;
                              TL_range = range(0.0, 50.0, length=200))
    ω_vec  = Float64[]
    Te_vec = Float64[]

    for TL in TL_range
        _, ia_ss, ω_ss, Te_ss = dc_ind_steady_state(p, Vf, Va, Float64(TL))
        push!(ω_vec,  ω_ss)
        push!(Te_vec, Te_ss)
    end

    return ω_vec, Te_vec
end

"""
    dc_ind_torque_speed_v2(p, Vf, Va; n_points) -> (ω_vec, Te_vec)

Genera curva par-velocidad barriendo velocidades desde 0 hasta velocidad en vacío.
Más preciso que barrer torques cuando la curva es casi horizontal.
"""
function dc_ind_torque_speed_v2(p::DCMotorIndParams, Vf::Float64, Va::Float64;
                                 n_points::Int = 300)
    KΦ   = p.Ke * (Vf / p.Rf)
    ω_0  = Va / KΦ                        # velocidad en vacío (TL=0, ia→0)
    ω_vec  = collect(range(0.0, ω_0 * 1.05, length=n_points))
    Te_vec = Float64[]

    for ω in ω_vec
        ia = (Va - KΦ * ω) / p.Ra
        Te = KΦ * ia
        push!(Te_vec, max(Te, 0.0))        # solo región motora
    end

    return ω_vec, Te_vec
end
