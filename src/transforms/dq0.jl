# =============================================================================
# Transformada de Park (dq0)
# Convierte magnitudes trifásicas (abc) a marco rotante (dq0)
# θ es el ángulo eléctrico del rotor en radianes
# =============================================================================

"""
    park(a, b, c, θ) -> (d, q, z)

Transformada de Park: convierte magnitudes trifásicas `a`, `b`, `c`
al marco de referencia rotante dq0 dado el ángulo eléctrico `θ` (rad).

# Ejemplo
```julia
vd, vq, v0 = park(va, vb, vc, θ)
```
"""
function park(a::T, b::T, c::T, θ::T) where T <: Real
    d =  (2/3) * ( a*cos(θ)       + b*cos(θ - 2π/3) + c*cos(θ + 2π/3))
    q = -(2/3) * ( a*sin(θ)       + b*sin(θ - 2π/3) + c*sin(θ + 2π/3))
    z =  (1/3) * ( a + b + c)
    return d, q, z
end

"""
    park_inverse(d, q, z, θ) -> (a, b, c)

Transformada de Park inversa: convierte componentes dq0
de vuelta a magnitudes trifásicas abc.
"""
function park_inverse(d::T, q::T, z::T, θ::T) where T <: Real
    a = d*cos(θ)             - q*sin(θ)             + z
    b = d*cos(θ - 2π/3)     - q*sin(θ - 2π/3)      + z
    c = d*cos(θ + 2π/3)     - q*sin(θ + 2π/3)      + z
    return a, b, c
end