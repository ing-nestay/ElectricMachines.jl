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
# Tests Motor DC
# =============================================================================
@testset "Motor DC" begin

    params = DCMotorParams(
        1.0,    # R = 1 Ω
        0.01,   # L = 10 mH
        0.1,    # Ke = 0.1 V·s/rad
        0.01,   # J = 0.01 kg·m²
        0.001   # B = 0.001 N·m·s/rad
    )

    V, TL = 12.0, 0.0

    # Estado estacionario analítico correcto
    ω_ss = (params.Ke*V - params.R*TL) / (params.Ke^2 + params.R*params.B)
    i_ss = (V - params.Ke*ω_ss) / params.R
    x_ss = [i_ss, ω_ss]

    # En SS las derivadas deben ser ≈ 0
    dx = [0.0, 0.0]
    dc_motor_ode!(dx, x_ss, params, 0.0, V, TL)

    @test abs(dx[1]) < 1e-8
    @test abs(dx[2]) < 1e-8

end
