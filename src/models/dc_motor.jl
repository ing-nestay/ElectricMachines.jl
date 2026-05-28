# =============================================================================
# Modelo Motor DC con excitación separada
# Estados: x = [i, ω]  (corriente de armadura, velocidad angular)
# Entradas: V (tensión), TL (torque de carga)
# =============================================================================

struct DCMotorParams
    R::Float64    # Resistencia de armadura [Ω]
    L::Float64    # Inductancia de armadura [H]
    Ke::Float64   # Constante de back-EMF / torque [V·s/rad]
    J::Float64    # Inercia del rotor [kg·m²]
    B::Float64    # Fricción viscosa [N·m·s/rad]
end

"""
    dc_motor_ode!(dx, x, params, t, V, TL)

Ecuaciones de estado del motor DC.
- x[1] = i  (corriente de armadura [A])
- x[2] = ω  (velocidad angular [rad/s])
"""
function dc_motor_ode!(dx, x, params::DCMotorParams, t, V, TL)
    i, ω = x[1], x[2]
    dx[1] = (V - params.R*i - params.Ke*ω) / params.L
    dx[2] = (params.Ke*i - params.B*ω - TL) / params.J
end

"""
    dc_steady_state(params, V, TL) -> (i_ss, ω_ss)

Calcula el punto de operación en estado estacionario.
"""
function dc_steady_state(params::DCMotorParams, V::Float64, TL::Float64)
    i_ss = (V - params.Ke * (V - params.R*(TL/params.Ke)) / params.R) / params.R
    ω_ss = (V - params.R*i_ss) / params.Ke
    return i_ss, ω_ss
end
