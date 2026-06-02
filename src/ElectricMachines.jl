module ElectricMachines

using LinearAlgebra
using StaticArrays

# Transformadas de referencia
include("transforms/dq0.jl")
include("transforms/alpha_beta.jl")
include("transforms/space_phasor.jl")

# Modelos de máquinas
include("models/inductances.jl")          # §1.3 inductancias propias y mutuas
include("models/park_equations.jl")       # §1.7 ecuaciones de Park
include("models/electromagnetic_torque.jl") # §1.8 momento electromagnético
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