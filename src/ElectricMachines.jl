module ElectricMachines

using LinearAlgebra
using StaticArrays

# Transformadas de referencia
include("transforms/dq0.jl")
include("transforms/alpha_beta.jl")

# Modelos de máquinas
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