"""
    visual_tests.jl

Pruebas visuales para verificar la correctitud de los módulos implementados.
Ejecutar con el entorno del proyecto activado:

    julia --project=. examples/visual_tests.jl

O en el REPL de VS Code (con Revise):

    include("examples/visual_tests.jl")
"""

using ElectricMachines
using Plots
using DifferentialEquations

# Motor de inducción de referencia (4 kW, 380 V, 50 Hz, 2 pares de polos)
# Parámetros de circuito equivalente típicos
motor = induction_motor_params_from_equivalent_circuit(
    Rs = 1.405,   # Ω
    Rr = 1.395,   # Ω
    Xs = 1.840,   # Ω  (dispersión estator)
    Xr = 1.840,   # Ω  (dispersión rotor)
    Xm = 40.0,    # Ω  (magnetización)
    J  = 0.0131,  # kg·m²
    p  = 2,
    fn = 50.0
)

ωs = motor.ωs  # 314.16 rad/s
Vs = 380.0 / sqrt(3)  # tensión de fase [V]

println("=" ^ 60)
println("Motor de referencia:")
println("  Rs=$(motor.Rs) Ω, Rr=$(motor.Rr) Ω")
println("  Ls=$(round(motor.Ls*1000,digits=2)) mH, Lr=$(round(motor.Lr*1000,digits=2)) mH")
println("  Lsr=$(round(motor.Lsr*1000,digits=2)) mH")
println("  σ = $(round(1 - motor.Lsr^2/(motor.Ls*motor.Lr), digits=4))")
println("=" ^ 60)

# =============================================================================
# FIGURA 1 — Curva par-velocidad (estado estacionario)
# Verifica: §1.9 isotropic_steady_state + §1.8 torque
# Criterio: forma de campana con máximo (par de arranque y par máximo visibles)
# =============================================================================
println("\n[1/4] Curva par-velocidad...")

# Convertir motor InductionMotorParams → IsotropicMachineParams
iso = IsotropicMachineParams(
    motor.Rs, motor.Rr, motor.Ls, motor.Lr, motor.Lsr,
    motor.J, motor.p, motor.ωs
)

slips    = range(0.001, 0.999, length=500)
torques  = Float64[]
speeds   = Float64[]

for s in slips
    _, _, _, _, Te = isotropic_steady_state(iso, complex(Vs), s)
    push!(torques, Te)
    push!(speeds,  ωs * (1 - s) / motor.p)  # velocidad mecánica [rad/s]
end

speeds_rpm = speeds .* 60 / (2π)

p1 = plot(speeds_rpm, torques,
    xlabel = "Velocidad [rpm]",
    ylabel = "Torque electromagnético [N·m]",
    title  = "Fig.1 — Curva Par-Velocidad (Estado Estacionario)",
    label  = "Te(n)",
    lw = 2, color = :blue,
    legend = :topleft)
hline!(p1, [0.0], color=:black, lw=0.5, label="")
vline!(p1, [3000.0], color=:gray, ls=:dash, lw=1, label="ωs")

# Marcar par máximo
idx_max = argmax(torques)
scatter!(p1, [speeds_rpm[idx_max]], [torques[idx_max]],
    color=:red, ms=6, label="Te_max = $(round(torques[idx_max],digits=1)) N·m")

println("  Par máximo: $(round(torques[idx_max],digits=2)) N·m a $(round(speeds_rpm[idx_max],digits=0)) rpm")
println("  Par de arranque (s=1): $(round(torques[end],digits=2)) N·m")

# =============================================================================
# FIGURA 2 — Arranque dinámico (transitorio desde ω=0)
# Verifica: induction_ode! + integración numérica
# Criterio: ω converge a ωs·(1-s_ss), Te converge a Tm en SS
# =============================================================================
println("\n[2/4] Arranque dinámico (transitorio)...")

Tm_load = 5.0   # N·m  (torque de carga constante)
vs_phasor = complex(Vs, 0.0)
vr_phasor = complex(0.0, 0.0)

# Condiciones iniciales: máquina en reposo, flujos en SS de baja velocidad
x0 = [0.0, -Vs/ωs, 0.0, 0.0, 0.0]  # ψs ≈ -jVs/ωs, ψr=0, ω=0

tspan = (0.0, 1.5)

prob = ODEProblem(
    (dx, x, p, t) -> induction_ode!(dx, x, motor, t, vs_phasor, vr_phasor, -Tm_load),
    x0, tspan
)
sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-8)

t_vec  = sol.t
ω_vec  = [sol.u[i][5] for i in eachindex(sol.u)]
ωm_vec = ω_vec ./ motor.p  # velocidad mecánica rad/s

# Reconstruir torque a lo largo de la solución
Te_vec = Float64[]
for i in eachindex(sol.u)
    x   = sol.u[i]
    ψs  = complex(x[1], x[2])
    is, _ = induction_currents_from_flux(ψs, complex(x[3], x[4]), motor)
    Te  = induction_torque(ψs, is, motor)
    push!(Te_vec, Te)
end

p2a = plot(t_vec, ωm_vec .* 60/(2π),
    xlabel="Tiempo [s]", ylabel="Velocidad [rpm]",
    title="Fig.2 — Arranque Dinámico",
    label="ωm(t)", lw=2, color=:blue)
hline!(p2a, [ωs/motor.p * 60/(2π)], color=:gray, ls=:dash, lw=1, label="ωs")

p2b = plot(t_vec, Te_vec,
    xlabel="Tiempo [s]", ylabel="Torque [N·m]",
    label="Te(t)", lw=2, color=:red)
hline!(p2b, [Tm_load], color=:black, ls=:dash, lw=1, label="Tm = $(Tm_load) N·m")

p2 = plot(p2a, p2b, layout=(2,1))

ω_final = ωm_vec[end] * 60/(2π)
println("  Velocidad final: $(round(ω_final, digits=1)) rpm  (ωs = $(round(ωs/motor.p*60/(2π),digits=1)) rpm)")
println("  Torque final: $(round(Te_vec[end], digits=2)) N·m  (Tm = $Tm_load N·m)")

# =============================================================================
# FIGURA 3 — Transformada de Park: ida y vuelta
# Verifica: abc_to_dq + dq_to_abc (§1.6)
# Criterio: señal reconstruida = señal original
# =============================================================================
println("\n[3/4] Transformada de Park (reversibilidad)...")

t_sig = range(0, 2/50, length=1000)  # 2 ciclos a 50 Hz
ω50   = 2π * 50.0
ia = @. sqrt(2) * 100.0 * cos(ω50 * t_sig)
ib = @. sqrt(2) * 100.0 * cos(ω50 * t_sig - 2π/3)
ic = @. sqrt(2) * 100.0 * cos(ω50 * t_sig + 2π/3)

ia_rec = similar(ia)
ib_rec = similar(ib)
ic_rec = similar(ic)

for (k, t) in enumerate(t_sig)
    γ = ω50 * t
    id, iq, i0 = abc_to_dq(ia[k], ib[k], ic[k], γ)
    a, b, c    = dq_to_abc(id, iq, i0, γ)
    ia_rec[k] = a; ib_rec[k] = b; ic_rec[k] = c
end

err_max = maximum(abs.(ia .- ia_rec))
println("  Error máximo de reconstrucción: $(round(err_max, sigdigits=3)) A  (debe ser < 1e-10)")

p3 = plot(t_sig.*1000, ia, label="ia original", lw=2, color=:blue,
    xlabel="Tiempo [ms]", ylabel="Corriente [A]",
    title="Fig.3 — Park: reconstrucción abc→dq→abc")
plot!(p3, t_sig.*1000, ia_rec, label="ia reconstruida", lw=1.5,
    color=:red, ls=:dash)

# =============================================================================
# FIGURA 4 — Corriente de estator en arranque (plano complejo)
# Verifica: fasor espacial durante el transitorio
# Criterio: espiral que converge al punto de operación nominal
# =============================================================================
println("\n[4/4] Trayectoria del fasor de corriente...")

is_traj = ComplexF64[]
for i in eachindex(sol.u)
    x  = sol.u[i]
    ψs = complex(x[1], x[2])
    ψr = complex(x[3], x[4])
    is, _ = induction_currents_from_flux(ψs, ψr, motor)
    push!(is_traj, is)
end

p4 = plot(real.(is_traj), imag.(is_traj),
    xlabel="Re{is} [A]", ylabel="Im{is} [A]",
    title="Fig.4 — Trayectoria fasor is durante arranque",
    label="is(t)", lw=1, color=:purple)
scatter!(p4, [real(is_traj[1])], [imag(is_traj[1])],
    color=:green, ms=6, label="t=0")
scatter!(p4, [real(is_traj[end])], [imag(is_traj[end])],
    color=:red, ms=6, label="SS")

# =============================================================================
# Mostrar todas las figuras
# =============================================================================
fig = plot(p1, p2, p3, p4, layout=(2,2), size=(1200, 900))
display(fig)
savefig(fig, "examples/visual_tests_output.png")

println("\n✓ Figuras guardadas en examples/visual_tests_output.png")
println("✓ Tests visuales completados.")
