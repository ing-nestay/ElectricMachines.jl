# =============================================================================
# Transformada de Clarke (αβ0)
# Convierte magnitudes trifásicas (abc) a coordenadas estacionarias (αβ0)
# =============================================================================

"""
    clarke(a, b, c) -> (α, β, γ)

Transformada de Clarke: convierte magnitudes trifásicas `a`, `b`, `c`
a componentes en marco estacionario αβ0.

# Ejemplo
```julia
α, β, γ = clarke(ia, ib, ic)
```
"""
function clarke(a::T, b::T, c::T) where T <: Real
    α = (2/3) * (a - 0.5*b - 0.5*c)
    β = (2/3) * (√3/2 * b - √3/2 * c)
    γ = (1/3) * (a + b + c)      # componente homopolar (cero)
    return α, β, γ
end

"""
    clarke_inverse(α, β, γ) -> (a, b, c)

Transformada de Clarke inversa: convierte componentes αβ0
de vuelta a magnitudes trifásicas abc.
"""
function clarke_inverse(α::T, β::T, γ::T) where T <: Real
    a = α + γ
    b = -0.5*α + (√3/2)*β + γ
    c = -0.5*α - (√3/2)*β + γ
    return a, b, c
end