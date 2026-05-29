# =============================================================================
# Modelo Dinámico Máquina de Inducción (MAS)
# Marco teórico: Mora Castro, A. — UTFSM, 2022
# Coordenadas: Sistema de Coordenadas Sincrónico (SCS)
# Variables de estado complejas (fasores espaciales)
#
# Estados: x = [ψsα, ψsβ, ψrα, ψrβ, ω]
#   ψs = ψsα + j·ψsβ  → enlace de flujo estator [Wb]
#   ψr = ψrα + j·ψrβ  → enlace de flujo rotor   [Wb]
#   ω                  → velocidad angular eléctrica [rad_elec/s]
#
# Parámetros: Rs, Rr, Ls, Lr, Lsr, J, p
# Entradas:   vs (tensión estator), vr (tensión rotor, =0 SCIM), Tm (torque mec.)
# =============================================================================

struct InductionMotorParams
    Rs::Float64    # Resistencia estator [Ω]
    Rr::Float64    # Resistencia rotor referida al estator [Ω]
    Ls::Float64    # Inductancia propia estator [H]
    Lr::Float64    # Inductancia propia rotor [H]
    Lsr::Float64   # Inductancia mutua estator-rotor [H]  (= Lrs para SCIM)
    J::Float64     # Inercia total [kg·m²]
    p::Int         # Número de pares de polos
    ωs::Float64    # Frecuencia sincrónica [rad/s]
end

"""
    induction_motor_params_from_nameplate(; Vn, fn, Pn, η, pf, p) -> InductionMotorParams

Estima parámetros aproximados a partir de datos de placa.
Útil para arranque rápido sin ensayos de caracterización.
"""
function induction_motor_params_from_nameplate(;
    Vn::Float64,   # Tensión de línea nominal [V]
    fn::Float64,   # Frecuencia nominal [Hz]
    Pn::Float64,   # Potencia nominal [W]
    η::Float64,    # Eficiencia nominal (0-1)
    pf::Float64,   # Factor de potencia nominal (0-1)
    p::Int         # Pares de polos
)
    ωs = 2π * fn
    Vs = Vn / sqrt(3)            # tensión de fase
    In = Pn / (sqrt(3) * Vn * pf * η)

    # Estimaciones empíricas (Krause / Chapman)
    Rs  = 0.04 * Vs / In
    Rr  = 0.04 * Vs / In
    Lm  = 0.95 * Vs / (ωs * In)  # inductancia de magnetización
    Ls  = Lm + 0.05 * Lm
    Lr  = Lm + 0.05 * Lm
    Lsr = Lm
    J   = 0.5 * (Pn / (ωs/p)^2)  # estimación inercia

    return InductionMotorParams(Rs, Rr, Ls, Lr, Lsr, J, p, ωs)
end

"""
    induction_coupling_factors(p::InductionMotorParams)

Retorna factores de acoplamiento kr, ks y coeficiente de dispersión σ
(ec. 1.33-1.34 del documento de referencia).
"""
function induction_coupling_factors(p::InductionMotorParams)
    kr = p.Lsr / p.Lr          # factor acoplamiento rotor
    ks = p.Lsr / p.Ls          # factor acoplamiento estator (Lrs=Lsr para SCIM)
    σ  = 1.0 - ks * kr         # coeficiente dispersión total
    return kr, ks, σ
end

"""
    induction_currents_from_flux(ψs, ψr, p) -> (is, ir)

Obtiene corrientes de estator y rotor a partir de los enlaces de flujo.
Invierte el sistema (ec. 1.31-1.32):
  ψs = Ls·is + Lsr·ir
  ψr = Lr·ir + Lsr·is
"""
function induction_currents_from_flux(ψs::ComplexF64, ψr::ComplexF64,
                                       p::InductionMotorParams)
    # Inversa analítica del sistema 2×2
    det = p.Ls * p.Lr - p.Lsr^2
    is  = (p.Lr * ψs - p.Lsr * ψr) / det
    ir  = (p.Ls * ψr - p.Lsr * ψs) / det
    return is, ir
end

"""
    induction_torque(ψs, is, p) -> Te

Momento electromagnético según ec. (1.42):
  Te = -(3/2)·p·Im{ψs · conj(is)}
"""
function induction_torque(ψs::ComplexF64, is::ComplexF64,
                           p::InductionMotorParams)
    return -(3/2) * p.p * imag(ψs * conj(is))
end

"""
    induction_ode!(dx, x, params, t, vs, vr, Tm)

Ecuaciones de estado de la MAS en coordenadas sincrónicas (SCS).

Estado x = [ψsα, ψsβ, ψrα, ψrβ, ω]:
- x[1:2] → Re e Im del fasor ψs (flujo estator)
- x[3:4] → Re e Im del fasor ψr (flujo rotor)
- x[5]   → velocidad angular eléctrica ω [rad_elec/s]

Entradas:
- vs::ComplexF64  tensión estator (fasor espacial)
- vr::ComplexF64  tensión rotor   (= 0 para SCIM jaula de ardilla)
- Tm::Float64     torque mecánico externo [N·m]
"""
function induction_ode!(dx, x, p::InductionMotorParams, t,
                         vs::ComplexF64, vr::ComplexF64, Tm::Float64)
    # Reconstruir fasores de flujo desde el vector de estado
    ψs = complex(x[1], x[2])
    ψr = complex(x[3], x[4])
    ω  = x[5]

    # Frecuencia de deslizamiento (ec. 1.39)
    ωr = p.ωs - ω

    # Corrientes desde flujos (ec. 1.31-1.32 invertidas)
    is, ir = induction_currents_from_flux(ψs, ψr, p)

    # Torque electromagnético (ec. 1.42)
    Te = induction_torque(ψs, is, p)

    # Ecuaciones diferenciales electromagnéticas (ec. 1.35-1.36)
    dψs = vs - p.Rs * is - im * p.ωs * ψs
    dψr = vr - p.Rr * ir - im * ωr  * ψr

    # Ecuación mecánica (ec. 1.55)
    dω  = (p.p / p.J) * (Te + Tm)

    # Escribir en el vector de derivadas
    dx[1] = real(dψs)
    dx[2] = imag(dψs)
    dx[3] = real(dψr)
    dx[4] = imag(dψr)
    dx[5] = dω
end

"""
    induction_steady_state(p, Vs_amp, Tm) -> (ψs_ss, ψr_ss, is_ss, ir_ss, ω_ss, Te_ss, s_ss)

Calcula el punto de operación en estado estacionario dado Vs y Tm.
Resuelve el sistema algebraico (derivadas = 0).
Retorna también el deslizamiento s = ωr/ωs.
"""
function induction_steady_state(p::InductionMotorParams,
                                 Vs_amp::Float64, Tm::Float64)
    # Tensión estator en eje d del SCS (orientada con vs)
    vs = complex(Vs_amp, 0.0)
    vr = complex(0.0, 0.0)   # SCIM: rotor cortocircuitado

    # Flujo de estator impuesto por la red (ec. 1.72): ψs = -j·Vs/ωs
    ψs_ss = -im * Vs_amp / p.ωs

    # Buscar ω_ss tal que Te = -Tm en SS (dω/dt = 0)
    # Barrer deslizamiento s ∈ (0, 1) y encontrar cruce
    best_ω = p.ωs * 0.95   # valor inicial
    best_err = Inf

    for s in range(0.001, 0.5, length=5000)
        ω_test = p.ωs * (1 - s)
        ωr_test = p.ωs * s

        # ψr desde ec. diferencial en SS (dψr/dt = 0):
        # 0 = vr - Rr·ir - j·ωr·ψr  →  ψr = (vr - Rr·ir) / (j·ωr)
        # Combinando con ec. de flujo: sistema lineal en ψr
        det  = p.Ls * p.Lr - p.Lsr^2
        A    = p.Rr * p.Ls / det + im * ωr_test
        B    = -p.Rr * p.Lsr / det
        ψr_t = (vr - B * ψs_ss) / A

        is_t, ir_t = induction_currents_from_flux(ψs_ss, ψr_t, p)
        Te_t = induction_torque(ψs_ss, is_t, p)

        err = abs(Te_t + Tm)
        if err < best_err
            best_err = err
            best_ω   = ω_test
            ψr_ss    = ψr_t
            is_ss    = is_t
            ir_ss    = ir_t
            Te_ss    = Te_t
        end
    end

    ψr_ss = complex(0.0)
    is_ss = complex(0.0)
    ir_ss = complex(0.0)
    Te_ss = 0.0

    for s in range(0.001, 0.5, length=5000)
        ω_test  = p.ωs * (1 - s)
        ωr_test = p.ωs * s
        det     = p.Ls * p.Lr - p.Lsr^2
        A       = p.Rr * p.Ls / det + im * ωr_test
        B       = -p.Rr * p.Lsr / det
        ψr_t    = (vr - B * ψs_ss) / A
        is_t, ir_t = induction_currents_from_flux(ψs_ss, ψr_t, p)
        Te_t    = induction_torque(ψs_ss, is_t, p)
        err     = abs(Te_t + Tm)
        if err < best_err
            best_err = err
            best_ω   = ω_test
            ψr_ss    = ψr_t
            is_ss    = is_t
            ir_ss    = ir_t
            Te_ss    = Te_t
        end
    end

    s_ss = (p.ωs - best_ω) / p.ωs
    return ψs_ss, ψr_ss, is_ss, ir_ss, best_ω, Te_ss, s_ss
end

# =============================================================================
# Versión corregida del estado estacionario
# =============================================================================
function induction_steady_state_v2(p::InductionMotorParams,
                                    Vs_amp::Float64, Tm::Float64)
    vs = complex(Vs_amp, 0.0)
    vr = complex(0.0, 0.0)

    # Flujo estator impuesto por la red (ec. 1.72)
    ψs_ss = -im * Vs_amp / p.ωs

    best_ω   = p.ωs * 0.95
    best_err = Inf
    ψr_ss    = complex(0.0)
    is_ss    = complex(0.0)
    ir_ss    = complex(0.0)
    Te_ss    = 0.0

    for s in range(0.0001, 0.5, length=10000)
        ω_test  = p.ωs * (1.0 - s)
        ωr_test = p.ωs * s

        # ψr en SS: resuelve 0 = vr - Rr·ir - j·ωr·ψr
        # con ir = (Ls·ψr - Lsr·ψs) / det
        det  = p.Ls * p.Lr - p.Lsr^2
        A    = p.Rr * p.Ls / det + im * ωr_test
        B    = -p.Rr * p.Lsr / det
        ψr_t = (vr - B * ψs_ss) / A

        is_t, ir_t = induction_currents_from_flux(ψs_ss, ψr_t, p)
        Te_t       = induction_torque(ψs_ss, is_t, p)

        err = abs(Te_t + Tm)
        if err < best_err
            best_err = err
            best_ω   = ω_test
            ψr_ss    = ψr_t
            is_ss    = is_t
            ir_ss    = ir_t
            Te_ss    = Te_t
        end
    end

    s_ss = (p.ωs - best_ω) / p.ωs
    return ψs_ss, ψr_ss, is_ss, ir_ss, best_ω, Te_ss, s_ss
end

# =============================================================================
# Constructor desde parámetros del circuito equivalente (ensayos)
# Notación estándar: Rs, Rr', Xs, Xr', Xm (todos en Ω a frecuencia nominal)
# Referencia: Chapman cap. 6, Krause cap. 2
# =============================================================================
"""
    induction_motor_params_from_equivalent_circuit(; Rs, Rr, Xs, Xr, Xm, J, p, fn)

Construye InductionMotorParams desde parámetros del circuito equivalente.
- Rs, Rr : resistencias estator y rotor [Ω]
- Xs, Xr : reactancias de dispersión estator y rotor [Ω]
- Xm     : reactancia de magnetización [Ω]
- J      : inercia [kg·m²]
- p      : pares de polos
- fn     : frecuencia nominal [Hz]
"""
function induction_motor_params_from_equivalent_circuit(;
    Rs::Float64, Rr::Float64,
    Xs::Float64, Xr::Float64, Xm::Float64,
    J::Float64,  p::Int, fn::Float64)

    ωs  = 2π * fn
    Lls = Xs / ωs    # inductancia dispersión estator
    Llr = Xr / ωs    # inductancia dispersión rotor
    Lm  = Xm / ωs    # inductancia de magnetización
    Ls  = Lls + Lm   # inductancia propia estator
    Lr  = Llr + Lm   # inductancia propia rotor
    Lsr = Lm         # inductancia mutua

    return InductionMotorParams(Rs, Rr, Ls, Lr, Lsr, J, p, ωs)
end
