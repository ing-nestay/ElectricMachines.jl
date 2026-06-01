# =============================================================================
# Inductancias propias y mutuas de máquinas de campo giratorio
#
# Derivadas desde primeros principios mediante:
#   - Ley de Ampere aplicada al entrehierro (campo unidimensional)
#   - Desarrollo en series de Fourier de la permeancia y la fmm
#   - Integración del flujo enlazado por cada devanado
#
# Notación:
#   γ   : ángulo eléctrico del rotor respecto al eje magnético de la fase a [rad]
#   L1  : valor medio de la inductancia propia de una fase del estator [H]
#   L2  : amplitud de variación de la inductancia propia (por anisotropía) [H]
#   Lad : inductancia propia máxima (eje d alineado con fase a)
#   Laq : inductancia propia mínima (eje q alineado con fase a)
#   L1f : valor máximo de la inductancia mutua estator-campo [H]
#   L1D : valor máximo de la inductancia mutua estator-amortiguador eje d [H]
#   L1Q : valor máximo de la inductancia mutua estator-amortiguador eje q [H]
# =============================================================================

# -----------------------------------------------------------------------------
# Inductancias propias del devanado de armadura (estator trifásico)
# Ec. (1.3.17)–(1.3.19): resultado de integrar la inducción producida por cada
# fase sobre la bobina concentrada equivalente de paso completo de la misma fase.
# La variación con γ tiene origen en la anisotropía magnética del rotor.
# -----------------------------------------------------------------------------

"""
    Laa(L1, L2, γ) -> Float64

Inductancia propia de la fase a.
"""
Laa(L1::Float64, L2::Float64, γ::Float64) = L1 + L2 * cos(2γ)

"""
    Lbb(L1, L2, γ) -> Float64

Inductancia propia de la fase b.
"""
Lbb(L1::Float64, L2::Float64, γ::Float64) = L1 + L2 * cos(2γ + 2π/3)

"""
    Lcc(L1, L2, γ) -> Float64

Inductancia propia de la fase c.
"""
Lcc(L1::Float64, L2::Float64, γ::Float64) = L1 + L2 * cos(2γ - 2π/3)

# -----------------------------------------------------------------------------
# Inductancias mutuas entre fases del devanado de armadura
# Ec. (1.3.21)–(1.3.24): resultado de integrar la inducción producida por una
# fase sobre la bobina equivalente de otra fase.
# El signo negativo refleja la disposición física a 120° entre ejes magnéticos.
# -----------------------------------------------------------------------------

"""
    Lab(L1, L2, γ) -> Float64

Inductancia mutua entre fases a y b.
"""
Lab(L1::Float64, L2::Float64, γ::Float64) = -L1/2 + L2 * cos(2γ - 2π/3)

"""
    Lac(L1, L2, γ) -> Float64

Inductancia mutua entre fases a y c.
"""
Lac(L1::Float64, L2::Float64, γ::Float64) = -L1/2 + L2 * cos(2γ + 2π/3)

"""
    Lbc(L1, L2, γ) -> Float64

Inductancia mutua entre fases b y c.
"""
Lbc(L1::Float64, L2::Float64, γ::Float64) = -L1/2 + L2 * cos(2γ)

# -----------------------------------------------------------------------------
# Inductancias mutuas estator–campo
# Ec. (1.3.28)–(1.3.31): el devanado de campo tiene su eje alineado con el eje d
# del rotor, por lo que la mutua es máxima cuando γ = 0 y nula cuando γ = π/2.
# -----------------------------------------------------------------------------

"""
    Lfa(L1f, γ) -> Float64

Inductancia mutua entre el devanado de campo y la fase a del estator.
"""
Lfa(L1f::Float64, γ::Float64) = L1f * cos(γ)

"""
    Lfb(L1f, γ) -> Float64

Inductancia mutua entre el devanado de campo y la fase b del estator.
"""
Lfb(L1f::Float64, γ::Float64) = L1f * cos(γ - 2π/3)

"""
    Lfc(L1f, γ) -> Float64

Inductancia mutua entre el devanado de campo y la fase c del estator.
"""
Lfc(L1f::Float64, γ::Float64) = L1f * cos(γ + 2π/3)

# -----------------------------------------------------------------------------
# Inductancias mutuas estator–amortiguador eje d (devanado D)
# Ec. (1.3.38): misma dependencia que el devanado de campo (eje d).
# -----------------------------------------------------------------------------

LDa(L1D::Float64, γ::Float64) = L1D * cos(γ)
LDb(L1D::Float64, γ::Float64) = L1D * cos(γ - 2π/3)
LDc(L1D::Float64, γ::Float64) = L1D * cos(γ + 2π/3)

# -----------------------------------------------------------------------------
# Inductancias mutuas estator–amortiguador eje q (devanado Q)
# Ec. (1.3.38): el eje q está en cuadratura con el eje d → dependencia senoidal.
# -----------------------------------------------------------------------------

LQa(L1Q::Float64, γ::Float64) = -L1Q * sin(γ)
LQb(L1Q::Float64, γ::Float64) = -L1Q * sin(γ - 2π/3)
LQc(L1Q::Float64, γ::Float64) = -L1Q * sin(γ + 2π/3)

# -----------------------------------------------------------------------------
# Matrices de inductancia completas
# Ec. (1.3.39)–(1.3.44): reúne todos los acoplamientos del sistema estator-rotor.
# El enlace de flujo de cada devanado = fila de la matriz × vector de corrientes.
# -----------------------------------------------------------------------------

"""
    inductance_matrix_stator(L1, L2, Lσ1, γ) -> Matrix{Float64}

Submatriz (3×3) de inductancias de la armadura (estator), incluyendo dispersión Lσ1.
Diagonal: inductancias propias; fuera de diagonal: mutuas entre fases.
"""
function inductance_matrix_stator(L1::Float64, L2::Float64,
                                   Lσ1::Float64, γ::Float64)
    return [Lσ1 + Laa(L1, L2, γ)   Lab(L1, L2, γ)          Lac(L1, L2, γ);
            Lab(L1, L2, γ)          Lσ1 + Lbb(L1, L2, γ)   Lbc(L1, L2, γ);
            Lac(L1, L2, γ)          Lbc(L1, L2, γ)          Lσ1 + Lcc(L1, L2, γ)]
end

"""
    inductance_matrix_stator_rotor(L1f, L1D, L1Q, γ) -> Matrix{Float64}

Submatriz (3×3) de inductancias mutuas estator–rotor [fases a,b,c] × [campo f, amort. D, amort. Q].
"""
function inductance_matrix_stator_rotor(L1f::Float64, L1D::Float64,
                                         L1Q::Float64, γ::Float64)
    return [Lfa(L1f, γ)   LDa(L1D, γ)   LQa(L1Q, γ);
            Lfb(L1f, γ)   LDb(L1D, γ)   LQb(L1Q, γ);
            Lfc(L1f, γ)   LDc(L1D, γ)   LQc(L1Q, γ)]
end

# -----------------------------------------------------------------------------
# Parámetros Lm, L1, L2 desde geometría del entrehierro
# Ec. (1.3.12)–(1.3.16): expresiones analíticas en función de dimensiones físicas.
# Útil para diseño o para obtener parámetros desde datos geométricos.
#
#   Lm   : inductancia de magnetización (rotor isotrópico, α = 1)
#   Lad  : inductancia eje d (máxima de la propia de armadura)
#   Laq  : inductancia eje q (mínima de la propia de armadura)
#   L1   : valor medio = (Lad + Laq) / 2
#   L2   : semivariación = (Lad - Laq) / 2
# -----------------------------------------------------------------------------

"""
    magnetizing_inductance(μ0, R, l, N1, fd1, δ_eff, p) -> Float64

Inductancia de magnetización Lm para rotor isotrópico (α = 1).
Parámetros: μ0 [H/m], R radio medio [m], l longitud axial [m],
N1 vueltas por fase, fd1 factor de distribución/cuerda, δ_eff entrehierro efectivo [m], p pares de polos.
"""
function magnetizing_inductance(μ0::Float64, R::Float64, l::Float64,
                                 N1::Float64, fd1::Float64,
                                 δ_eff::Float64, p::Int)
    return (4μ0 / π^2) * (R * l / δ_eff) * (N1 * fd1 / p)^2
end

"""
    salient_pole_inductances(Lm, α) -> (Lad, Laq, L1, L2)

Calcula inductancias Lad, Laq, L1, L2 para un rotor de polos salientes
dado el factor de arco polar α (relación ancho zapata / paso polar, 0 < α ≤ 1).
α = 1 corresponde al rotor cilíndrico isotrópico (Lad = Laq = Lm).
"""
function salient_pole_inductances(Lm::Float64, α::Float64)
    Lad = Lm * (α + sin(α * π) / π)
    Laq = Lm * (α - sin(α * π) / π)
    L1  = (Lad + Laq) / 2
    L2  = (Lad - Laq) / 2
    return Lad, Laq, L1, L2
end
