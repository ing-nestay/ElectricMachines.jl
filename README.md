# ElectricMachines.jl

[![Build Status](https://github.com/ing-nestay/ElectricMachines.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ing-nestay/ElectricMachines.jl/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Julia](https://img.shields.io/badge/Julia-1.10+-purple.svg)](https://julialang.org)

> Librería Julia para modelado y simulación de máquinas eléctricas de campo giratorio.  
> Formulación desde primeros principios físicos: Ampere, Faraday, Lorentz y transformación de Park.

---

## Motivación

La mayoría de las herramientas de simulación de máquinas eléctricas son propietarias (MATLAB/Simulink).  
Esta librería ofrece una alternativa **abierta, reproducible y extensible** para investigadores,  
estudiantes e ingenieros, aprovechando el rendimiento numérico de Julia.

Toda la formulación se deriva explícitamente desde primeros principios, sin dependencia de fuentes bibliográficas específicas:

- Ley de Ampere aplicada al entrehierro (campo unidimensional)
- Series de Fourier de la permeancia y la fmm
- Ley de Faraday para los enlaces de flujo
- Fasor espacial como variable compleja
- Transformación de Park (coordenadas fijas al rotor)
- Fuerzas de Lorentz para el momento electromagnético

---

## Modelos implementados

| Módulo | Archivo | Estado |
|--------|---------|--------|
| Fasor espacial (§1.6) | `src/transforms/space_phasor.jl` | ✅ |
| Transformada de Clarke (αβ0) | `src/transforms/alpha_beta.jl` | ✅ |
| Transformada de Park (dq0) | `src/transforms/dq0.jl` | ✅ |
| Inductancias propias y mutuas (§1.3) | `src/models/inductances.jl` | ✅ |
| Ecuaciones de Park (§1.7) | `src/models/park_equations.jl` | ✅ |
| Torque electromagnético (§1.8) | `src/models/electromagnetic_torque.jl` | ✅ |
| Máquina isotrópica simétrica (§1.9) | `src/models/isotropic_machine.jl` | ✅ |
| Motor de inducción trifásico (MAS) | `src/models/induction.jl` | ✅ |
| Motor DC — Excitación independiente | `src/models/dc_motor.jl` | ✅ |
| Motor síncrono (MS) | `src/models/synchronous.jl` | 🔄 En desarrollo |

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
using ElectricMachines

# Clarke: abc → αβ0
α, β, γ = clarke(1.0, -0.5, -0.5)

# Park: abc → dq0
vd, vq, v0 = park(1.0, -0.5, -0.5, π/4)

# Fasor espacial de corriente
is = spatial_phasor(ia, ib, ic)   # is = (2/3)(ia + a·ib + a²·ic)
```

### Motor de inducción trifásico (MAS)

```julia
using ElectricMachines
using DifferentialEquations

# Parámetros desde circuito equivalente (ensayos de vacío y cortocircuito)
motor = induction_motor_params_from_equivalent_circuit(
    Rs = 1.405,   # Ω
    Rr = 1.395,   # Ω
    Xs = 1.840,   # Ω
    Xr = 1.840,   # Ω
    Xm = 40.0,    # Ω
    J  = 0.0131,  # kg·m²
    p  = 2,
    fn = 50.0
)

# Estado estacionario para un deslizamiento dado
Vs  = 380.0 / sqrt(3)
iso = IsotropicMachineParams(motor.Rs, motor.Rr, motor.Ls, motor.Lr,
                              motor.Lsr, motor.J, motor.p, motor.ωs)
is, ir, ψs, ψr, Te = isotropic_steady_state(iso, complex(Vs), 0.05)

# Arranque dinámico (simulación ODE)
vs = complex(Vs, 0.0)
x0 = [0.0, -Vs/motor.ωs, 0.0, 0.0, 0.0]
prob = ODEProblem(
    (dx, x, p, t) -> induction_ode!(dx, x, motor, t, vs, 0.0+0.0im, -5.0),
    x0, (0.0, 1.5)
)
sol = solve(prob, Tsit5())
```

### Inductancias de máquina síncrona con polos salientes

```julia
using ElectricMachines

# Inductancias desde geometría del rotor
Lm = 100e-3       # inductancia de magnetización [H]
α  = 0.7          # factor de arco polar (relación zapata/paso polar)
Lad, Laq, L1, L2 = salient_pole_inductances(Lm, α)

# Matriz de inductancias del estator para posición γ del rotor
M = inductance_matrix_stator(L1, L2, 10e-3, γ)
```

---

## Base matemática

### Ecuaciones de estado de la MAS (coordenadas sincrónicas)

Las variables de estado son los flujos de enlace $\psi_s$ y $\psi_r$ (fasores espaciales):

$$\frac{d\psi_s}{dt} = v_s - R_s \, i_s - j\omega_s \, \psi_s$$

$$\frac{d\psi_r}{dt} = v_r - R_r \, i_r - j\omega_r \, \psi_r \quad (\omega_r = \omega_s - \omega)$$

$$\frac{d\omega}{dt} = \frac{p}{J}\left(T_e - T_m\right)$$

El torque electromagnético por fuerzas de Lorentz (§1.8):

$$T_e = -\frac{3}{2} p \, \text{Im}\left\{\psi_s \cdot \bar{i}_s\right\}$$

### Inductancias propias y mutuas (§1.3)

Las inductancias dependen del ángulo eléctrico $\gamma$ del rotor:

$$L_{aa}(\gamma) = L_1 + L_2 \cos(2\gamma)$$

$$L_{ab}(\gamma) = -\frac{L_1}{2} + L_2 \cos\!\left(2\gamma - \frac{2\pi}{3}\right)$$

Para rotor isotrópico ($L_2 = 0$): $L_{aa} = L_{bb} = L_{cc} = L_s$ (constante).

---

## Ejemplos y tests visuales

```julia
# Tests visuales de la MAS (curva par-vel., arranque, Park, fasor)
include("examples/visual_tests.jl")

# Tests visuales del módulo de inductancias
include("examples/visual_tests_inductances.jl")
```

| Script | Figuras generadas |
|--------|-------------------|
| `visual_tests.jl` | Curva par-velocidad, arranque dinámico, reversibilidad de Park, trayectoria fasor |
| `visual_tests_inductances.jl` | Inductancias propias, mutuas, estator-rotor, positividad definida, saliencia |

---

## Tests unitarios

```bash
julia --project=test test/runtests.jl
```

```
Test Summary:                       | Pass  Total
Transformada Clarke                 |    6      6
Transformada Park                   |    6      6
Motor DC — Excitación Independiente |    6      6
```

---

## Convenciones

| Símbolo | Descripción |
|---------|-------------|
| $\psi$ | Enlace de flujo [Wb] — variable de estado principal |
| $\omega$ | Velocidad eléctrica [rad_elec/s] |
| $\omega_m = \omega/p$ | Velocidad mecánica [rad/s] |
| $\gamma$ | Ángulo eléctrico rotor–eje magnético fase a [rad] |
| $p$ | Número de pares de polos |
| $s = \omega_r/\omega_s$ | Deslizamiento |
| $T > 0$ | Torque motor (acelera), $T_m > 0$ carga resistente |

La transformada de Park usa el factor **2/3** (no normalizado).  
La componente de secuencia cero: $i_0 = (i_a + i_b + i_c)/3$, nula en estrella sin neutro.

---

## Estructura del repositorio

```
src/
├── ElectricMachines.jl          # módulo principal (exports)
├── transforms/
│   ├── space_phasor.jl          # §1.6 fasor espacial
│   ├── alpha_beta.jl            # transformada de Clarke
│   └── dq0.jl                   # transformada de Park
├── models/
│   ├── inductances.jl           # §1.3 inductancias
│   ├── park_equations.jl        # §1.7 ecuaciones de Park
│   ├── electromagnetic_torque.jl # §1.8 torque
│   ├── isotropic_machine.jl     # §1.9 máquina isotrópica
│   ├── induction.jl             # MAS dinámica
│   ├── dc_motor.jl              # motor DC
│   └── synchronous.jl           # MS (en desarrollo)
examples/
├── visual_tests.jl              # tests visuales MAS
└── visual_tests_inductances.jl  # tests visuales inductancias
test/
└── runtests.jl                  # tests unitarios
```

---

## Licencia

MIT © [ing-nestay](https://github.com/ing-nestay)
