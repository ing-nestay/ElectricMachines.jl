using Test

# Ruta base del proyecto (sube un nivel desde /test)
const PROJECT_ROOT = joinpath(@__DIR__, "..")

# Cargar transformadas
include(joinpath(PROJECT_ROOT, "src/transforms/alpha_beta.jl"))
include(joinpath(PROJECT_ROOT, "src/transforms/dq0.jl"))

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