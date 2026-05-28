using Plots
using DifferentialEquations

include("../src/models/dc_motor.jl")

p = DCMotorIndParams(
    150.0, 10.0, 1.0, 0.01, 0.8, 0.05, 0.01
)

Vf       = 220.0
k_load   = 0.0008   # TL = k·ω²  [N·m·s²/rad²]
colores  = [:dodgerblue, :darkorange, :green, :purple]
tensiones = [220.0, 180.0, 140.0, 100.0]

# =============================================================================
# FIGURA 1: Curvas par-velocidad + carga cuadrática + puntos de operación
# =============================================================================
plt1 = plot(
    title   = "Curvas Par-Velocidad — Motor DC Excitación Independiente",
    xlabel  = "Velocidad ω [rad/s]",
    ylabel  = "Par [N·m]",
    legend  = :topright,
    grid    = true,
    size    = (860, 540),
    xlims   = (0, 320),
    ylims   = (0, 70)
)

puntos_op = []

for (Va_i, col) in zip(tensiones, colores)
    ω_vec, Te_vec = dc_ind_torque_speed_v2(p, Vf, Va_i)
    plot!(plt1, ω_vec, Te_vec,
          label     = "Motor Va = $(Int(Va_i)) V",
          linewidth = 2.5,
          color     = col)

    # Intersección analítica: k·ω² + (KΦ²/Ra + B)·ω - KΦ·Va/Ra = 0
    KΦ    = p.Ke * (Vf / p.Rf)
    a_c   = k_load
    b_c   = KΦ^2 / p.Ra + p.B
    c_c   = -KΦ * Va_i / p.Ra
    disc  = b_c^2 - 4*a_c*c_c
    if disc >= 0
        ω_op = (-b_c + sqrt(disc)) / (2*a_c)
        if ω_op > 0
            Te_op = k_load * ω_op^2
            push!(puntos_op, (ω_op, Te_op, Va_i, col))
        end
    end
end

# Curva de carga cuadrática
ω_load = range(0.0, 310.0, length=400)
plot!(plt1, collect(ω_load), k_load .* collect(ω_load).^2,
      label     = "Carga: TL = k·ω²",
      linewidth = 2,
      linestyle = :dash,
      color     = :black)

# Puntos de operación
for (ω_op, Te_op, Va_i, col) in puntos_op
    scatter!(plt1, [ω_op], [Te_op],
             markersize  = 9,
             markercolor = col,
             markerstrokecolor = :black,
             markerstrokewidth = 1,
             label = "OP $(Int(Va_i))V → ω=$(round(ω_op,digits=1)) rad/s, Te=$(round(Te_op,digits=1)) N·m")
end

display(plt1)
savefig(plt1, "examples/torque_speed_quadratic_load.png")
println("✓ Figura 1 guardada")

# =============================================================================
# FIGURA 2: Arranque dinámico con carga cuadrática — 4 paneles
# =============================================================================
x0     = [0.0, 0.0, 0.0]
t_span = (0.0, 8.0)
Va     = 220.0

ode_quad!(dx, x, p, t) = dc_ind_ode!(dx, x, p, t, Vf, Va, k_load * x[3]^2)
prob = ODEProblem(ode_quad!, x0, t_span, p)
sol  = solve(prob, Tsit5(), saveat=0.02)

t_  = sol.t
i_f = [u[1] for u in sol.u]
i_a = [u[2] for u in sol.u]
ω_  = [u[3] for u in sol.u]
Te_ = [p.Ke * u[1] * u[2] for u in sol.u]
TL_ = k_load .* ω_ .^ 2

plt2 = plot(layout=(4,1), size=(820, 920),
            grid=true, linewidth=2, left_margin=8Plots.mm, bottom_margin=3Plots.mm)

plot!(plt2[1], t_, i_f, color=:dodgerblue, label=false,
      title="Corriente de campo  if [A]", ylabel="if [A]")

plot!(plt2[2], t_, i_a, color=:red, label=false,
      title="Corriente de armadura  ia [A]", ylabel="ia [A]")

plot!(plt2[3], t_, ω_, color=:green, label=false,
      title="Velocidad angular  ω [rad/s]", ylabel="ω [rad/s]")

plot!(plt2[4], t_, Te_, label="Te (motor)", color=:purple, linewidth=2)
plot!(plt2[4], t_, TL_, label="TL = k·ω²", color=:darkorange,
      linestyle=:dash, linewidth=2,
      title="Par electromagnético vs carga  [N·m]",
      ylabel="T [N·m]", xlabel="Tiempo [s]")

display(plt2)
savefig(plt2, "examples/startup_quadratic_load.png")
println("✓ Figura 2 guardada")
