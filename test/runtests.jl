using Test

const PROJECT_ROOT = joinpath(@__DIR__, "..")

# Cargar módulos
include(joinpath(PROJECT_ROOT, "src/transforms/alpha_beta.jl"))
include(joinpath(PROJECT_ROOT, "src/transforms/dq0.jl"))
include(joinpath(PROJECT_ROOT, "src/models/dc_motor.jl"))

# =============================================================================
# Tests Transformada de Clarke (αβ)
# =============================================================================
@testset "Transformada Clarke" begin

    α, β, γ = clarke(1.0, -0.5, -0.5)
    @test α ≈ 1.0  atol=1e-10
    @test β ≈ 0.0  atol=1e-10
    @test γ ≈ 0.0  atol=1e-10

    a, b, c = 1.0, -0.5, -0.5
    α, β, γ = clarke(a, b, c)
    a2, b2, c2 = clarke_inverse(α, β, γ)
    @test a2 ≈ a  atol=1e-10
    @test b2 ≈ b  atol=1e-10
    @test c2 ≈ c  atol=1e-10

end

# =============================================================================
# Tests Transformada de Park (dq0)
# =============================================================================
@testset "Transformada Park" begin

    vd, vq, v0 = park(1.0, -0.5, -0.5, 0.0)
    @test vd ≈ 1.0  atol=1e-10
    @test vq ≈ 0.0  atol=1e-10
    @test v0 ≈ 0.0  atol=1e-10

    a, b, c, θ = 1.0, -0.5, -0.5, π/4
    d, q, z = park(a, b, c, θ)
    a2, b2, c2 = park_inverse(d, q, z, θ)
    @test a2 ≈ a  atol=1e-10
    @test b2 ≈ b  atol=1e-10
    @test c2 ≈ c  atol=1e-10

end

# =============================================================================
# Tests Motor DC — Excitación Independiente
# =============================================================================
@testset "Motor DC — Excitación Independiente" begin

    p = DCMotorIndParams(
        150.0,   # Rf [Ω]
        10.0,    # Lf [H]
        1.0,     # Ra [Ω]
        0.01,    # La [H]
        0.8,     # Ke [V·s/A·rad]
        0.05,    # J  [kg·m²]
        0.01     # B  [N·m·s/rad]
    )

    Vf, Va, TL = 220.0, 220.0, 5.0

    # Test 1: derivadas ≈ 0 en estado estacionario
    if_ss, ia_ss, ω_ss, Te_ss = dc_ind_steady_state(p, Vf, Va, TL)
    dx = zeros(3)
    dc_ind_ode!(dx, [if_ss, ia_ss, ω_ss], p, 0.0, Vf, Va, TL)
    @test abs(dx[1]) < 1e-8
    @test abs(dx[2]) < 1e-8
    @test abs(dx[3]) < 1e-8

    # Test 2: balance de potencia eléctrica = mecánica + pérdidas
    P_elec = Va * ia_ss
    P_mec  = Te_ss * ω_ss
    P_Ra   = p.Ra * ia_ss^2
    @test P_elec ≈ P_mec + P_Ra  atol=1e-6

    # Test 3: curva par-velocidad tiene pendiente negativa
    ω_vec, Te_vec = dc_ind_torque_speed(p, Vf, Va, TL_range=range(0.0, 20.0, length=50))
    @test length(ω_vec) == 50
    @test ω_vec[1] > ω_vec[end]

end