# =============================================================================
# Máquina isotrópica simétrica — caso particular (§1.9)
#
# Derivación:
#   Cuando el entrehierro es constante (rotor cilíndrico, α = 1) desaparece
#   la anisotropía magnética: L2 = 0, Lad = Laq = Lm.
#   Las inductancias propias y mutuas del estator se vuelven constantes
#   (independientes de γ) y los ejes d y q son equivalentes:
#       L1d = L1q = Lσ1 + (3/2)·Lm  ≡  Ls
#
#   Las ecuaciones de Park se simplifican a:
#       vs,r = R1·is,r + (j·ω + d/dt)·ψs,r          (ec. 1.7.9, L2=0)
#
#   En coordenadas sincrónicas (ω = ωs = cte) el sistema se hace lineal.
#   Esta es la base del modelo dinámico de la MAS (jaula de ardilla) y de
#   la MS de rotor cilíndrico.
#
#   Para la MAS: vr = 0, ωr = ωs − ω (frecuencia de deslizamiento)
#   Para la MS cilíndrica: vr = vf (campo de excitación), ωr = 0 en SS
# =============================================================================

"""
    IsotropicMachineParams

Parámetros del modelo isotrópico simétrico (rotor cilíndrico).
Válido para MAS y MS de rotor liso.
"""
struct IsotropicMachineParams
    R1::Float64    # resistencia estator [Ω]
    Rr::Float64    # resistencia rotor referida al estator [Ω]
    Ls::Float64    # inductancia propia estator: Lσ1 + (3/2)·Lm [H]
    Lr::Float64    # inductancia propia rotor referida [H]
    Lm::Float64    # inductancia de magnetización: (3/2)·L1f [H]
    J::Float64     # inercia total [kg·m²]
    p::Int         # pares de polos
    ωs::Float64    # frecuencia sincrónica [rad_elec/s]
end

# -----------------------------------------------------------------------------
# Flujos de enlace (sistema isotrópico, ec. 1.7.11 con L1d = L1q = Ls)
#   ψs = Ls·is + Lm·ir
#   ψr = Lr·ir + Lm·is
# Válido como ecuación compleja en cualquier sistema de referencia girante.
# -----------------------------------------------------------------------------

"""
    isotropic_flux_linkages(is, ir, p) -> (ψs, ψr)

Enlace de flujo para la máquina isotrópica en variables complejas.
    ψs = Ls·is + Lm·ir
    ψr = Lr·ir + Lm·is
"""
function isotropic_flux_linkages(is::ComplexF64, ir::ComplexF64,
                                  p::IsotropicMachineParams)
    ψs = p.Ls * is + p.Lm * ir
    ψr = p.Lr * ir + p.Lm * is
    return ψs, ψr
end

"""
    isotropic_currents_from_flux(ψs, ψr, p) -> (is, ir)

Invierte el sistema de flujos de la máquina isotrópica.
    det = Ls·Lr − Lm²
    is  = (Lr·ψs − Lm·ψr) / det
    ir  = (Ls·ψr − Lm·ψs) / det
"""
function isotropic_currents_from_flux(ψs::ComplexF64, ψr::ComplexF64,
                                       p::IsotropicMachineParams)
    det = p.Ls * p.Lr - p.Lm^2
    is  = (p.Lr * ψs - p.Lm * ψr) / det
    ir  = (p.Ls * ψr - p.Lm * ψs) / det
    return is, ir
end

# -----------------------------------------------------------------------------
# Ecuaciones de estado en coordenadas sincrónicas (ω_ref = ωs)
# Para MAS: vr = 0, frecuencia de deslizamiento ωr = ωs − ω
# La linealidad para ω = cte permite soluciones analíticas vía Laplace.
#
#   dψs/dt = vs − R1·is − j·ωs·ψs
#   dψr/dt = vr − Rr·ir − j·ωr·ψr      con ωr = ωs − ω
#   dω/dt  = (p/J)·(T − Tm)
# -----------------------------------------------------------------------------

"""
    isotropic_ode!(dx, x, p, t, vs, vr, Tm)

Ecuaciones de estado de la máquina isotrópica en coordenadas sincrónicas.

Estado x = [Re(ψs), Im(ψs), Re(ψr), Im(ψr), ω]:
- x[1:2] → fasor ψs (flujo estator)
- x[3:4] → fasor ψr (flujo rotor)
- x[5]   → velocidad eléctrica ω [rad_elec/s]

Para MAS jaula: vr = 0+0j.
Para MS cilíndrica: vr = tensión de campo referida.
"""
function isotropic_ode!(dx, x, p::IsotropicMachineParams, t,
                         vs::ComplexF64, vr::ComplexF64, Tm::Float64)
    ψs = complex(x[1], x[2])
    ψr = complex(x[3], x[4])
    ω  = x[5]

    ωr = p.ωs - ω

    is, ir = isotropic_currents_from_flux(ψs, ψr, p)
    Te     = torque_phasor(ψs, is, p.p)

    dψs = vs - p.R1 * is - im * p.ωs * ψs
    dψr = vr - p.Rr * ir - im * ωr   * ψr
    dω  = (p.p / p.J) * (Te - Tm)

    dx[1] = real(dψs);  dx[2] = imag(dψs)
    dx[3] = real(dψr);  dx[4] = imag(dψr)
    dx[5] = dω
end

# -----------------------------------------------------------------------------
# Estado estacionario simétrico (§1.9, velocidad constante)
# En SS con velocidad constante dω/dt = 0 y las derivadas de flujo son
# puramente rotacionales:  dψ/dt → j·ω_ref·ψ
#
# Sistema algebraico complejo resultante:
#   vs = (R1 + j·ωs·Ls)·is + j·ωs·Lm·ir
#    0 = (Rr + j·ωr·Lr)·ir + j·ωr·Lm·is
# Que corresponde al circuito equivalente monofásico clásico.
# -----------------------------------------------------------------------------

"""
    isotropic_steady_state(p, Vs, s) -> (is, ir, ψs, ψr, Te)

Calcula el estado estacionario para un deslizamiento dado s = ωr/ωs.
Resuelve el sistema algebraico complejo (derivadas = 0 en SS).

Entradas:
- Vs : tensión de fase del estator (fasor, típicamente real) [V]
- s  : deslizamiento (0 = sincrónico, 1 = rotor bloqueado)
"""
function isotropic_steady_state(p::IsotropicMachineParams,
                                 Vs::ComplexF64, s::Float64)
    ωr = s * p.ωs

    # Sistema 2×2 complejo: [vs; 0] = Z · [is; ir]
    Z11 = p.R1 + im * p.ωs * p.Ls
    Z12 = im * p.ωs * p.Lm
    Z21 = im * ωr   * p.Lm
    Z22 = p.Rr + im * ωr   * p.Lr

    # Solución por regla de Cramer
    det_Z = Z11 * Z22 - Z12 * Z21
    is    = (Vs * Z22) / det_Z
    ir    = (-Vs * Z21) / det_Z

    ψs, ψr = isotropic_flux_linkages(is, ir, p)
    Te     = torque_phasor(ψs, is, p.p)

    return is, ir, ψs, ψr, Te
end

# -----------------------------------------------------------------------------
# Factor de dispersión total (coeficiente de Blondel)
# σ = 1 − Lm²/(Ls·Lr)
# Determina la inductancia transitoria y la dinámica rápida.
# -----------------------------------------------------------------------------

"""
    leakage_factor(p::IsotropicMachineParams) -> Float64

Coeficiente de dispersión total σ = 1 − Lm²/(Ls·Lr).
Relacionado con los factores de acoplamiento: σ = 1 − ks·kr.
"""
leakage_factor(p::IsotropicMachineParams) = 1.0 - p.Lm^2 / (p.Ls * p.Lr)
