# ElectricMachines.jl

[![Build Status](https://github.com/ing-nestay/ElectricMachines.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ing-nestay/ElectricMachines.jl/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Julia](https://img.shields.io/badge/Julia-1.9+-purple.svg)](https://julialang.org)

> Librería Julia para simulación dinámica de máquinas eléctricas rotativas.  
> Modelos en espacio de estados, transformadas de referencia y análisis de régimen permanente.

---

## Motivación

La mayoría de las herramientas de simulación de máquinas eléctricas son propietarias (MATLAB/Simulink).  
Esta librería busca ofrecer una alternativa **abierta, reproducible y extensible** para investigadores,  
estudiantes e ingenieros, aprovechando el rendimiento numérico de Julia.

---

## Modelos implementados

| Modelo | Archivo | Estado |
|--------|---------|--------|
| Transformada de Clarke (αβ) | `src/transforms/alpha_beta.jl` | ✅ |
| Transformada de Park (dq0) | `src/transforms/dq0.jl` | ✅ |
| Motor DC — Excitación Independiente | `src/models/dc_motor.jl` | ✅ |
| Motor DC — Excitación Shunt | `src/models/dc_motor.jl` | 🔄 En desarrollo |
| Motor DC — Excitación Compound | `src/models/dc_motor.jl` | 📋 Planificado |
| Motor de Inducción Trifásico | `src/models/induction.jl` | 📋 Planificado |
| Motor Síncrono PMSM | `src/models/synchronous.jl` | 📋 Planificado |

---

## Instalación

```julia
using Pkg
Pkg.develop(url="https://github.com/ing-nestay/ElectricMachines.jl")
```

---

## Uso rápido

### Transformadas de referencia

```julia
include("src/transforms/alpha_beta.jl")
include("src/transforms/dq0.jl")

# Clarke: trifásico → αβ
α, β, γ = clarke(1.0, -0.5, -0.5)

# Park: trifásico → dq0 con ángulo θ
vd, vq, v0 = park(1.0, -0.5, -0.5, π/4)
```

### Motor DC — Excitación Independiente

```julia
include("src/models/dc_motor.jl")

# Definir parámetros
p = DCMotorIndParams(
    150.0,  # Rf [Ω]
    10.0,   # Lf [H]
    1.0,    # Ra [Ω]
    0.01,   # La [H]
    0.8,    # Ke [V·s/A·rad]
    0.05,   # J  [kg·m²]
    0.01    # B  [N·m·s/rad]
)

# Punto de operación en régimen permanente
if_ss, ia_ss, ω_ss, Te_ss = dc_ind_steady_state(p, 220.0, 220.0, 5.0)

# Curva par-velocidad
ω_vec, Te_vec = dc_ind_torque_speed_v2(p, 220.0, 220.0)
```

---

## Modelo matemático

### Motor DC Excitación Independiente

Tres ecuaciones diferenciales acopladas:

$$L_f \frac{di_f}{dt} = V_f - R_f \, i_f$$

$$L_a \frac{di_a}{dt} = V_a - R_a \, i_a - K_e \, i_f \, \omega$$

$$J \frac{d\omega}{dt} = K_e \, i_f \, i_a - B \, \omega - T_L$$

El par electromagnético es $T_e = K_e \cdot i_f \cdot i_a$.

---

## Ejemplos

| Ejemplo | Descripción |
|---------|-------------|
| [`examples/dc_independent_excitation.jl`](examples/dc_independent_excitation.jl) | Curvas par-velocidad y arranque dinámico con carga cuadrática |

---

## Tests

```bash
julia --project=. test/runtests.jl
```Test Summary:                       | Pass  Total
Transformada Clarke                 |    6      6
Transformada Park                   |    6      6
Motor DC — Excitación Independiente |    6      6

---

## Referencias

- Chapman, S. *Electric Machinery Fundamentals*, 5th ed. McGraw-Hill, 2012.
- Krause, P.C. *Analysis of Electric Machinery and Drive Systems*, 3rd ed. IEEE Press, 2013.
- Mohan, N. *Electric Drives*, 3rd ed. Wiley, 2012.

---

## Contribuir

¿Quieres agregar un modelo? Lee [`CONTRIBUTING.md`](CONTRIBUTING.md) *(próximamente)*.  
Abre un [issue](https://github.com/ing-nestay/ElectricMachines.jl/issues) para reportar errores o proponer mejoras.

---

## Licencia

MIT © [ing-nestay](https://github.com/ing-nestay)