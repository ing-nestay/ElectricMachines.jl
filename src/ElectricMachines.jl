module ElectricMachines

using LinearAlgebra
using StaticArrays

# --- Transformadas ---
export clarke, clarke_inverse
export park, park_inverse
export spatial_phasor, phasor_to_rotating, phasor_from_rotating
export zero_sequence, phasor_to_abc, phasor_to_alphabeta, alphabeta_to_phasor
export abc_to_dq, dq_to_abc

# --- Inductancias ---
export Laa, Lbb, Lcc, Lab, Lac, Lbc
export Lfa, Lfb, Lfc, LDa, LDb, LDc, LQa, LQb, LQc
export inductance_matrix_stator, inductance_matrix_stator_rotor
export magnetizing_inductance, salient_pole_inductances

# --- Ecuaciones de Park ---
export ParkInductances
export park_inductances, park_flux_linkages, park_flux_derivatives, park_currents_from_flux

# --- Torque electromagnético ---
export torque_dq, torque_phasor, mechanical_ode, electromagnetic_power

# --- Máquina isotrópica ---
export IsotropicMachineParams
export isotropic_flux_linkages, isotropic_currents_from_flux
export isotropic_ode!, isotropic_steady_state, leakage_factor

# --- Motor de inducción ---
export InductionMotorParams
export induction_ode!, induction_steady_state, induction_steady_state_v2
export induction_currents_from_flux, induction_torque
export induction_motor_params_from_nameplate, induction_motor_params_from_equivalent_circuit

# --- Motor DC ---
export DCMotorIndParams, dc_ind_ode!, dc_ind_steady_state, dc_ind_torque_speed

# Transformadas de referencia
include("transforms/dq0.jl")
include("transforms/alpha_beta.jl")
include("transforms/space_phasor.jl")

# Modelos de máquinas
include("models/inductances.jl")          # §1.3 inductancias propias y mutuas
include("models/park_equations.jl")       # §1.7 ecuaciones de Park
include("models/electromagnetic_torque.jl") # §1.8 momento electromagnético
include("models/isotropic_machine.jl")     # §1.9 máquina isotrópica simétrica
include("models/dc_motor.jl")
include("models/induction.jl")
include("models/synchronous.jl")

# Solvers
include("solvers/state_space.jl")

# Análisis
include("analysis/torque_speed.jl")

# Parámetros
include("parameters/machine_data.jl")

end # module ElectricMachines