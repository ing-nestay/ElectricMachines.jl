"""
    visual_tests_inductances.jl

Pruebas visuales para el módulo `models/inductances.jl`.
Verifica propiedades físicas esperadas desde primeros principios:
  - Periodicidad y simetría de inductancias propias y mutuas
  - Simetría de la matriz de inductancias (reciprocidad)
  - Positividad definida de la matriz estator
  - Comportamiento límite: rotor isotrópico (L2=0) y polos salientes
  - Inductancia de magnetización desde geometría

Ejecutar en REPL con entorno activado:
    include("examples/visual_tests_inductances.jl")
"""

using ElectricMachines
using Plots
using LinearAlgebra

println("=" ^ 60)
println("Tests visuales — Módulo de Inductancias (§1.3)")
println("=" ^ 60)

# Parámetros de referencia (máquina síncrona con polos salientes)
L1  = 100e-3   # [H]  valor medio inductancia propia
L2  = 30e-3    # [H]  semivariación (anisotropía)
Lσ1 = 10e-3    # [H]  dispersión estator
L1f = 80e-3    # [H]  mutua estator-campo
L1D = 75e-3    # [H]  mutua estator-amortiguador D
L1Q = 60e-3    # [H]  mutua estator-amortiguador Q

γ_vec = range(0, 2π, length=1000)  # barrido de posición del rotor

# =============================================================================
# FIGURA 1 — Inductancias propias de las tres fases vs ángulo del rotor
# Criterio físico:
#   - Período 2π en γ (pero variación con cos(2γ) → período π)
#   - Las tres curvas son idénticas pero desplazadas 2π/3 en argumento de cos(2γ)
#   - L2=0 → curvas constantes (rotor isotrópico)
#   - Siempre positivas (Laa > 0 para cualquier γ)
# =============================================================================
println("\n[1/5] Inductancias propias de las tres fases...")

Laa_v = [Laa(L1, L2, γ) for γ in γ_vec]
Lbb_v = [Lbb(L1, L2, γ) for γ in γ_vec]
Lcc_v = [Lcc(L1, L2, γ) for γ in γ_vec]

@assert all(Laa_v .> 0) "Laa debe ser siempre positiva"
@assert all(Lbb_v .> 0) "Lbb debe ser siempre positiva"
@assert all(Lcc_v .> 0) "Lcc debe ser siempre positiva"

# Verificar período π (la variación es cos(2γ)) — evaluando en puntos exactos
@assert isapprox(Laa(L1, L2, 0.0),        Laa(L1, L2, Float64(π)),     atol=1e-12) "Laa debe tener período π"
@assert isapprox(Laa(L1, L2, Float64(π)/4), Laa(L1, L2, 5Float64(π)/4), atol=1e-12) "Laa debe tener período π"

println("  Laa ∈ [$(round(minimum(Laa_v)*1e3,digits=1)), $(round(maximum(Laa_v)*1e3,digits=1))] mH")
println("  Período verificado: Laa(0) ≈ Laa(π) ✓")
println("  Todas positivas ✓")

p1 = plot(γ_vec, Laa_v.*1e3, label="Laa", lw=2, color=:blue,
    xlabel="Ángulo rotórico γ [rad]",
    ylabel="Inductancia [mH]",
    title="Fig.1 — Inductancias propias vs γ")
plot!(p1, γ_vec, Lbb_v.*1e3, label="Lbb", lw=2, color=:red)
plot!(p1, γ_vec, Lcc_v.*1e3, label="Lcc", lw=2, color=:green)
hline!(p1, [L1*1e3], color=:black, ls=:dash, lw=1, label="L1 (valor medio)")
vline!(p1, [π/2, π, 3π/2], color=:gray, ls=:dot, lw=0.8, label="")
xticks!(p1, ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3π/2", "2π"]))

# =============================================================================
# FIGURA 2 — Inductancias mutuas entre fases vs ángulo del rotor
# Criterio físico:
#   - Siempre negativas para L2=0 (fases a 120° → acoplamiento negativo)
#   - Con L2>0 pueden volverse positivas en ciertos ángulos
#   - Las tres curvas desplazadas 2π/3 entre sí
#   - Valor medio = -L1/2
# =============================================================================
println("\n[2/5] Inductancias mutuas entre fases...")

Lab_v = [ElectricMachines.Lab(L1, L2, γ) for γ in γ_vec]
Lac_v = [ElectricMachines.Lac(L1, L2, γ) for γ in γ_vec]
Lbc_v = [ElectricMachines.Lbc(L1, L2, γ) for γ in γ_vec]

# Valor medio debe ser -L1/2
mean_Lab = sum(Lab_v) / length(Lab_v)
@assert isapprox(mean_Lab, -L1/2, atol=1e-4) "Valor medio de Lab debe ser -L1/2"

println("  Valor medio Lab: $(round(mean_Lab*1e3,digits=2)) mH  (esperado: $(-L1/2*1e3) mH) ✓")
println("  Lab ∈ [$(round(minimum(Lab_v)*1e3,digits=1)), $(round(maximum(Lab_v)*1e3,digits=1))] mH")

p2 = plot(γ_vec, Lab_v.*1e3, label="Lab", lw=2, color=:blue,
    xlabel="Ángulo rotórico γ [rad]",
    ylabel="Inductancia [mH]",
    title="Fig.2 — Inductancias mutuas estator vs γ")
plot!(p2, γ_vec, Lac_v.*1e3, label="Lac", lw=2, color=:red)
plot!(p2, γ_vec, Lbc_v.*1e3, label="Lbc", lw=2, color=:green)
hline!(p2, [-L1/2*1e3], color=:black, ls=:dash, lw=1, label="-L1/2 (valor medio)")
hline!(p2, [0.0], color=:gray, lw=0.5, label="")
xticks!(p2, ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3π/2", "2π"]))

# =============================================================================
# FIGURA 3 — Inductancias mutuas estator–rotor (campo y amortiguadores)
# Criterio físico:
#   - Lfa: coseno, máxima en γ=0 (eje d alineado con fase a)
#   - LQa: seno negativo (eje q en cuadratura con eje d)
#   - Lfb, Lfc desplazadas ±2π/3 respecto a Lfa
#   - Lfa(0) = L1f (valor máximo de la mutua)
# =============================================================================
println("\n[3/5] Inductancias mutuas estator-rotor...")

Lfa_v = [Lfa(L1f, γ) for γ in γ_vec]
Lfb_v = [Lfb(L1f, γ) for γ in γ_vec]
Lfc_v = [Lfc(L1f, γ) for γ in γ_vec]
LQa_v = [LQa(L1Q, γ) for γ in γ_vec]

@assert isapprox(Lfa(L1f, 0.0), L1f,  atol=1e-12) "Lfa(0) debe ser L1f"
@assert isapprox(Lfa(L1f, π/2), 0.0,  atol=1e-12) "Lfa(π/2) debe ser 0"
@assert isapprox(LQa(L1Q, 0.0), 0.0,  atol=1e-12) "LQa(0) debe ser 0 (eje q ⊥ fase a en γ=0)"
@assert isapprox(LQa(L1Q, π/2), -L1Q, atol=1e-12) "LQa(π/2) debe ser -L1Q"

println("  Lfa(0) = $(round(Lfa(L1f,0.0)*1e3,digits=1)) mH = L1f ✓")
println("  Lfa(π/2) = $(round(Lfa(L1f,π/2)*1e3,digits=1)) mH = 0 ✓")
println("  LQa(0) = $(round(LQa(L1Q,0.0)*1e3,digits=1)) mH = 0 (cuadratura) ✓")

p3 = plot(γ_vec, Lfa_v.*1e3, label="Lfa (campo, fase a)", lw=2, color=:blue,
    xlabel="Ángulo rotórico γ [rad]",
    ylabel="Inductancia [mH]",
    title="Fig.3 — Inductancias mutuas estator–rotor vs γ")
plot!(p3, γ_vec, Lfb_v.*1e3, label="Lfb (campo, fase b)", lw=2, color=:red, ls=:dash)
plot!(p3, γ_vec, Lfc_v.*1e3, label="Lfc (campo, fase c)", lw=2, color=:green, ls=:dash)
plot!(p3, γ_vec, LQa_v.*1e3, label="LQa (amort. Q, fase a)", lw=2, color=:purple)
hline!(p3, [0.0], color=:gray, lw=0.5, label="")
xticks!(p3, ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3π/2", "2π"]))

# =============================================================================
# FIGURA 4 — Matriz de inductancias: simetría y positividad definida
# Criterio físico:
#   - Matriz simétrica para todo γ (reciprocidad electromagnética)
#   - Todos los autovalores positivos para todo γ (energía magnética > 0)
# =============================================================================
println("\n[4/5] Simetría y positividad definida de la matriz estator...")

γ_test = range(0, 2π, length=200)
eigvals_min = Float64[]
asym_error  = Float64[]

for γ in γ_test
    M = inductance_matrix_stator(L1, L2, Lσ1, γ)
    push!(eigvals_min, minimum(eigvals(M)))
    push!(asym_error,  maximum(abs.(M .- M')))
end

@assert all(eigvals_min .> 0) "La matriz debe ser definida positiva para todo γ"
@assert all(asym_error .< 1e-15) "La matriz debe ser simétrica para todo γ"

println("  Autovalor mínimo: $(round(minimum(eigvals_min)*1e3, digits=3)) mH > 0 ✓")
println("  Error de simetría máximo: $(maximum(asym_error)) ✓")

p4a = plot(γ_test, eigvals_min.*1e3,
    xlabel="Ángulo rotórico γ [rad]",
    ylabel="λ_min [mH]",
    title="Fig.4 — Autovalores mínimos de L_estator(γ)",
    label="λ_min(γ)", lw=2, color=:blue)
hline!(p4a, [0.0], color=:red, ls=:dash, lw=1, label="límite positividad")
xticks!(p4a, ([0, π/2, π, 3π/2, 2π], ["0", "π/2", "π", "3π/2", "2π"]))

# =============================================================================
# FIGURA 5 — Efecto de la saliencia: Lad, Laq, L1, L2 vs factor de arco polar α
# Criterio físico:
#   - α=1 (rotor cilíndrico): Lad = Laq = Lm, L2 = 0
#   - α→0 (saliencia extrema): Laq → 0, Lad → Lm·(1/π)·sin(π) → 0 también
#   - Lad > Laq siempre para 0 < α < 1
#   - L2 = (Lad-Laq)/2 > 0 siempre para α < 1
# =============================================================================
println("\n[5/5] Inductancias de polos salientes vs factor de arco polar α...")

Lm = 100e-3   # inductancia de magnetización isotrópica [H]
α_vec = range(0.1, 1.0, length=300)

Lad_v = Float64[]; Laq_v = Float64[]
L1_v  = Float64[]; L2_v  = Float64[]

for α in α_vec
    lad, laq, l1, l2 = salient_pole_inductances(Lm, α)
    push!(Lad_v, lad); push!(Laq_v, laq)
    push!(L1_v, l1);   push!(L2_v, l2)
end

# Verificar límite α=1 (rotor cilíndrico)
Lad1, Laq1, L1_1, L2_1 = salient_pole_inductances(Lm, 1.0)
@assert isapprox(Lad1, Laq1, atol=1e-10) "Para α=1: Lad debe ser igual a Laq"
@assert isapprox(L2_1, 0.0,  atol=1e-10) "Para α=1: L2 debe ser cero"
@assert all(Lad_v .>= Laq_v)             "Lad debe ser ≥ Laq para todo α"

println("  α=1.0 → Lad=$(round(Lad1*1e3,digits=2)) mH, Laq=$(round(Laq1*1e3,digits=2)) mH, L2=$(round(L2_1*1e3,digits=2)) mH ✓")
println("  Lad ≥ Laq para todo α ✓")

p5 = plot(α_vec, Lad_v.*1e3, label="Lad (eje d)", lw=2, color=:blue,
    xlabel="Factor de arco polar α",
    ylabel="Inductancia [mH]",
    title="Fig.5 — Saliencia: inductancias vs α")
plot!(p5, α_vec, Laq_v.*1e3, label="Laq (eje q)", lw=2, color=:red)
plot!(p5, α_vec, L1_v.*1e3,  label="L1 = (Lad+Laq)/2", lw=2, color=:green, ls=:dash)
plot!(p5, α_vec, L2_v.*1e3,  label="L2 = (Lad-Laq)/2", lw=2, color=:purple, ls=:dash)
hline!(p5, [Lm*1e3], color=:black, ls=:dot, lw=1, label="Lm (isotrópico)")
vline!(p5, [1.0], color=:gray, ls=:dot, lw=1, label="α=1 (cilíndrico)")

# =============================================================================
# Componer y mostrar figura final
# =============================================================================
fig = plot(p1, p2, p3, p4a, p5,
    layout = @layout([a b; c d; e{0.5w}]),
    size   = (1200, 1100))

display(fig)
savefig(fig, "examples/visual_tests_inductances_output.png")

println("\n✓ Figura guardada en examples/visual_tests_inductances_output.png")
println("✓ Todos los tests de inductancias completados.")
