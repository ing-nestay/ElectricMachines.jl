using Plots
using DifferentialEquations

include("../src/models/dc_motor.jl")

# =============================================================================
# Parámetros motor DC 1 kW excitación independiente
# =============================================================================
p = DCMotorIndParams(
    150.0,   # Rf [Ω]
    10.0,    # Lf [H]
    1.0,     # Ra [Ω]
    0.01,    # La [H]
    0.8,     # Ke [V·s/A·rad]
    0.05,    # J  [kg·m²]
    0.01     # B  [N·m·s/rad]
)

Vf, Va = 220.0, 220.0

# Coeficiente de carga cuadrática: TL = k_load * ω²
# Calibrado para que la intersección quede dentro del rango de operación
k_load = 0.0008   # [N·m·s²/rad²]

# =============================================================================
# FIGURA 1: Curva Par-Velocidad + curva de carga cuadrática
# =============================================================================
TL_range = range(0.0, 15.0, length=300)

plt1 = plot(
    title  = "Curvas Par-Velocidad — Motor DC Excitación Independiente",
    xlabel = "Velocidad ω [rad/s]",
    ylabel = "Par [N·m]",
    legend = :topright,
    grid   = true,
    size   = (800, 520)
)

# Curvas del motor para distintas tensiones
intersections = []
for Va_i in [220.0, 180.0, 140.0, 100.0]
    ω_vec, Te_vec = dc_ind_torque_speed(p, Vf, Va_i, TL_range=TL_range)
    plot!(plt1, ω_vec, Te_vec,
          label = "Motor Va = $(Int(Va_i)) V",
          linewidth = 2)

    # Calcular intersección: Te(ω) = k_load·ω²
    # en SS: Te_ss = k_load·ω_ss² y también Te_ss = KΦ·(Va - KΦ·ω_ss)/Ra
    KΦ = p.Ke * (Vf / p.Rf)
    # KΦ·Va/Ra - (KΦ²/Ra)·ω = k_load·ω²
    # k_load·ω² + (KΦ²/Ra)·ω - KΦ·Va_i/Ra = 0
    a_coef = k_load
    b_coef = KΦ^2 / p.Ra + p.B
    c_coef = -KΦ * Va_i / p.Ra
    discriminant = b_coef^2 - 4*a_coef*c_coef
    if discriminant >= 0
        ω_int = (-b_coef + sqrt(discriminant)) / (2*a_coef)
        Te_int = k_load * ω_int^2
        if ω_int > 0
            push!(intersections, (ω_int, Te_int, Va_i))
        end
    end
end

# Curva de carga cuadrática TL = k·ω²
ω_load = range(0.0, 300.0, length=300)
TL_load = k_load .* ω_load .^ 2
plot!(plt1, collect(ω_load), TL_load,
      label     = "Carga: TL = k·ω²",
      linewidth = 2,
      linestyle = :dash,
      color     = :black)

# Marcar puntos de intersección
for (ω_i, Te_i, Va_i) in intersections
    scatter!(plt1, [ω_i], [Te_i],
             markersize  = 8,
             markercolor = :red,
             markershape = :circle,
             label       = "OP Va=$(Int(Va_i)): ω=$(round(ω_i,digits=1)) rad/s, Te=$(round(Te_i,digits=2)) N·m")
end

display(plt1)
savefig(plt1, "examples/torque_speed_quadratic_load.png")
println("Figura 1 guardada")

# =============================================================================
# FIGURA 2: Arranque dinámico con carga cuadrática
# =============================================================================
x0     = [0.0, 0.0, 0.0]
t_span = (0.0, 8.0)

# Carga cuadrática depende del estado: TL(t) = k_load * ω(t)²
ode_quad!(dx, x, p, t) = dc_ind_ode!(dx, x, p, t, Vf, Va, k_load * x[3]^2)

prob = ODEProblem(ode_quad!, x0, t_span, p)
sol  = solve(prob, Tsit5(), saveat=0.01)

t_  = sol.t
i_f = [u[1] for u in sol.u]
i_a = [u[2] for u in sol.u]
ω_  = [u[3] for u in sol.u]
Te_ = [p.Ke * u[1] * u[2] for u in sol.u]
TL_ = k_load .* ω_ .^ 2

plt2 = plot(layout=(4,1), size=(800, 900), grid=true, linewidth=2,
            left_margin=5Plots.mm)

plot!(plt2[1], t_, i_f,
      title="Corriente de campo  if [A]", label=false, color=:blue)

plot!(plt2[2], t_, i_a,
      title="Corriente de armadura  ia [A]", label=false, color=:red)

plot!(plt2[3], t_, ω_,
      title="Velocidad angular  ω [rad/s]", label=false, color=:green)

plot!(plt2[4], t_, Te_, label="Te — motor", color=:purple, linewidth=2)
plot!(plt2[4], t_, TL_, label="TL — carga k·ω²",
      color=:orange, linestyle=:dash, linewidth=2,
      title="Par: electromagnético vs carga [N·m]",
      xlabel="Tiempo [s]")

display(plt2)
savefig(plt2, "examples/startup_quadratic_load.png")
println("Figura 2 guardada")
