# =============================================================================
# Tests Motor DC Excitación Independiente
# =============================================================================
include(joinpath(PROJECT_ROOT, "src/models/dc_motor.jl"))

@testset "Motor DC — Excitación Independiente" begin

    # Parámetros motor DC típico 1 kW
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

    # --- Test 1: estado estacionario analítico ---
    if_ss, ia_ss, ω_ss, Te_ss = dc_ind_steady_state(p, Vf, Va, TL)

    # Verificar con ODE: derivadas ≈ 0 en SS
    dx = zeros(3)
    dc_ind_ode!(dx, [if_ss, ia_ss, ω_ss], p, 0.0, Vf, Va, TL)
    @test abs(dx[1]) < 1e-8   # dif/dt ≈ 0
    @test abs(dx[2]) < 1e-8   # dia/dt ≈ 0
    @test abs(dx[3]) < 1e-8   # dω/dt  ≈ 0

    # --- Test 2: balance de potencia ---
    # Potencia eléctrica entrada ≈ Potencia mecánica + pérdidas
    P_elec = Va * ia_ss
    P_mec  = Te_ss * ω_ss
    P_Ra   = p.Ra * ia_ss^2
    @test P_elec ≈ P_mec + P_Ra  atol=1e-6

    # --- Test 3: curva par-velocidad tiene pendiente negativa ---
    ω_vec, Te_vec = dc_ind_torque_speed(p, Vf, Va, TL_range=range(0.0, 20.0, length=50))
    @test length(ω_vec) == 50
    @test ω_vec[1] > ω_vec[end]   # a mayor torque, menor velocidad

end
