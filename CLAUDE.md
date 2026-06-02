# ElectricMachines.jl — Contexto para Claude Code

## Propósito
Librería Julia para modelado y simulación de máquinas eléctricas de campo giratorio
(MAS y MS). Toda la formulación se deriva desde primeros principios físicos para
evitar dependencia de fuentes bibliográficas específicas.

## Base teórica
El fundamento analítico proviene del desarrollo estándar de máquinas de campo
giratorio basado en:
- Ley de Ampere aplicada al entrehierro (campo unidimensional)
- Desarrollo en series de Fourier de la permeancia y la fmm
- Ley de Faraday para los enlaces de flujo
- Fasor espacial como variable compleja (componente simétrica instantánea)
- Transformación de Park (coordenadas fijas al rotor)
- Fuerzas de Lorentz para el momento electromagnético

## Estructura del repositorio

```
src/
├── ElectricMachines.jl          # módulo principal
├── transforms/
│   ├── dq0.jl                   # transformada de Park (abc → dq0)
│   ├── alpha_beta.jl            # transformada de Clarke (abc → αβ0)
│   └── space_phasor.jl          # §1.6 fasor espacial completo
├── models/
│   ├── inductances.jl           # §1.3 inductancias propias y mutuas
│   ├── park_equations.jl        # §1.7 ecuaciones de Park (coordenadas dq)
│   ├── electromagnetic_torque.jl # §1.8 torque por fuerzas de Lorentz
│   ├── isotropic_machine.jl     # §1.9 máquina isotrópica simétrica
│   ├── induction.jl             # modelo dinámico MAS (coordenadas sincrónicas)
│   ├── induction_v2.jl          # respaldo versión anterior de induction.jl
│   ├── synchronous.jl           # modelo MS (vacío, pendiente)
│   ├── synchronous_v2.jl        # respaldo
│   └── dc_motor.jl              # motor DC excitación independiente
├── solvers/
│   └── state_space.jl
├── analysis/
│   └── torque_speed.jl
└── parameters/
    └── machine_data.jl
```

## Módulos implementados

### `transforms/space_phasor.jl` (§1.6)
Fasor espacial de corriente: `is = (2/3)(ia + a·ib + a²·ic)` con `a = e^(j2π/3)`.

Funciones principales:
- `spatial_phasor(ia, ib, ic)` → `ComplexF64`
- `phasor_to_rotating(is, γ)` → `is · e^(-jγ)`
- `phasor_from_rotating(is_r, γ)`
- `zero_sequence(ia, ib, ic)`
- `phasor_to_abc(is, i0)`
- `phasor_to_alphabeta(is)` / `alphabeta_to_phasor(iα, iβ)`
- `abc_to_dq(ia, ib, ic, γ)` / `dq_to_abc(id, iq, i0, γ)`

### `models/inductances.jl` (§1.3)
Inductancias propias y mutuas desde la integral de la inducción en el entrehierro.

Funciones principales:
- `Laa/Lbb/Lcc(L1, L2, γ)` — propias del estator
- `Lab/Lac/Lbc(L1, L2, γ)` — mutuas entre fases
- `Lfa/Lfb/Lfc(L1f, γ)` — mutuas estator–campo
- `LDa/LDb/LDc`, `LQa/LQb/LQc` — mutuas estator–amortiguadores
- `inductance_matrix_stator(L1, L2, Lσ1, γ)` → Matrix 3×3
- `inductance_matrix_stator_rotor(L1f, L1D, L1Q, γ)` → Matrix 3×3
- `magnetizing_inductance(μ0, R, l, N1, fd1, δ_eff, p)` — Lm desde geometría
- `salient_pole_inductances(Lm, α)` → `(Lad, Laq, L1, L2)`

### `models/park_equations.jl` (§1.7)
Ecuaciones de Park en forma de estado (flujos como variables de estado).

Tipos:
- `ParkInductances` — agrupa L1d, L1q, Lf, LD, LQ, Lf1, LD1, LQ1, LfD

Funciones principales:
- `park_inductances(...)` — calcula ParkInductances desde parámetros de entrehierro
- `park_flux_linkages(i1d, i1q, if_, iD, iQ, li)` → `(ψ1d, ψ1q, ψf, ψD, ψQ)`
- `park_flux_derivatives(v1d, v1q, vf, vD, vQ, ...)` → `(dψ1d, dψ1q, dψf, dψD, dψQ)`
- `park_currents_from_flux(ψ1d, ψ1q, ψf, ψD, ψQ, li)` → `(i1d, i1q, if_, iD, iQ)`

### `models/electromagnetic_torque.jl` (§1.8)
Momento electromagnético desde la integral de Lorentz.

Funciones principales:
- `torque_dq(ψ1d, ψ1q, i1d, i1q, p)` → `T = (3/2)·p·(ψ1d·i1q − ψ1q·i1d)`
- `torque_phasor(ψs, is, p)` → `T = −(3/2)·p·Im{ψs·conj(is)}`
- `mechanical_ode(T, Tm, ω, J, p)` → `(dω, dγ)`
- `electromagnetic_power(T, ω, p)`

### `models/isotropic_machine.jl` (§1.9)
Caso particular: rotor cilíndrico (L2 = 0, L1d = L1q = Ls).

Tipos:
- `IsotropicMachineParams` — R1, Rr, Ls, Lr, Lm, J, p, ωs

Funciones principales:
- `isotropic_flux_linkages(is, ir, p)`
- `isotropic_currents_from_flux(ψs, ψr, p)`
- `isotropic_ode!(dx, x, p, t, vs, vr, Tm)` — coordenadas sincrónicas
- `isotropic_steady_state(p, Vs, s)` — solución algebraica compleja
- `leakage_factor(p)` → `σ = 1 − Lm²/(Ls·Lr)`

### `models/induction.jl`
Modelo dinámico MAS en coordenadas sincrónicas (variables complejas).
Encabezado referencia Faraday (ec. 1.4.1), fasor espacial (ec. 1.6.4),
Park (§1.7) y Lorentz (§1.8) — sin cita a autores específicos.

Tipos: `InductionMotorParams` (Rs, Rr, Ls, Lr, Lsr, J, p, ωs)

Funciones: `induction_ode!`, `induction_steady_state`, `induction_steady_state_v2`,
`induction_currents_from_flux`, `induction_torque`,
`induction_motor_params_from_nameplate`, `induction_motor_params_from_equivalent_circuit`

## Convenciones

- Variables de estado: flujos de enlace ψ [Wb] (no corrientes)
- Velocidad: ω en rad_elec/s (eléctrica), ωm = ω/p en rad/s (mecánica)
- Ángulo γ: eléctrico, entre eje d del rotor y eje magnético de la fase a
- Factor 2/3 (no normalizado) en la transformación de Clarke/Park
- Componente de secuencia cero i0 = (ia+ib+ic)/3; nula para estrella sin neutro
- Torque motor positivo (T > 0 acelera), carga resistente Tm > 0

## Pendiente
- `synchronous.jl`: modelo completo MS (rotor saliente y cilíndrico)
- §1.10: excitación asimétrica y componentes simétricas
- Circuitos equivalentes en estado estacionario (MAS y MS)
- Tests unitarios para todos los módulos

## Rama de desarrollo activa
`claude/optimistic-dirac-UGeO9`

## Notas de sesión
- El token de GitHub para push debe configurarse como variable de entorno
  `GITHUB_TOKEN` en el entorno de Claude Code (Settings → Environments)
  para evitar pegarlo manualmente en el chat.
- Los archivos `*_v2.jl` son respaldos de versiones anteriores, no se incluyen
  en el módulo principal.
