# =============================================================================
# Ecuaciones de Park — coordenadas fijas al rotor (§1.7)
#
# Derivación:
#   Las tres ecuaciones de equilibrio del estator (ec. 1.4.2–1.4.4) se reducen
#   a una sola ecuación compleja mediante el fasor espacial (§1.6):
#
#       vs = R1·is + dψs/dt                          (ec. 1.7.1, marco estator)
#
#   Al transformar al sistema fijo al rotor (is = is,r · e^(jγ)), el operador
#   diferencial d/dt se convierte en (j·dγ/dt + d/dt), lo que introduce
#   la tensión de velocidad y linealiza el sistema para velocidad constante:
#
#       vs,r = R1·is,r + (j·ω + d/dt)·ψs,r          (ec. 1.7.9)
#
#   Proyectando sobre los ejes d y q (ec. 1.7.13–1.7.14):
#
#       v1d = R1·i1d + L1d·di1d/dt + L1D·diD/dt + L1f·dif/dt − ω·(L1q·i1q + L1Q·iQ)
#       v1q = R1·i1q + L1q·di1q/dt + L1Q·diQ/dt     + ω·(L1d·i1d + L1D·iD + L1f·if)
#
#   Ecuaciones del rotor (ec. 1.7.17–1.7.19):
#       vf  = Rf·if  + Lf·dif/dt  + Lf1·di1d/dt + LfD·diD/dt
#       vD  = RD·iD  + LD·diD/dt  + LD1·di1d/dt + LDf·dif/dt
#       vQ  = RQ·iQ  + LQ·diQ/dt  + LQ1·di1q/dt
#
#   Inductancias de campo giratorio (ec. 1.7.7 y 1.7.20):
#       LG1 = (3/2)·L1,  LG2 = (3/2)·L2
#       L1d = Lσ1 + LG1 + LG2,  L1q = Lσ1 + LG1 − LG2
#       Lf1 = (3/2)·L1f,  LD1 = (3/2)·L1D,  LQ1 = (3/2)·L1Q
# =============================================================================

"""
    ParkInductances

Inductancias del modelo de dos ejes (Park).
Agrupa las inductancias propias de los ejes d y q y las mutuas
estator–rotor referidas al marco girante.
"""
struct ParkInductances
    L1d::Float64   # inductancia propia eje d del estator [H]
    L1q::Float64   # inductancia propia eje q del estator [H]
    Lf::Float64    # inductancia propia del campo [H]
    LD::Float64    # inductancia propia amortiguador eje d [H]
    LQ::Float64    # inductancia propia amortiguador eje q [H]
    Lf1::Float64   # inductancia mutua campo–estator eje d [H]  = (3/2)·L1f
    LD1::Float64   # inductancia mutua amort.D–estator eje d [H] = (3/2)·L1D
    LQ1::Float64   # inductancia mutua amort.Q–estator eje q [H] = (3/2)·L1Q
    LfD::Float64   # inductancia mutua campo–amortiguador d [H]
end

"""
    park_inductances(L1, L2, Lσ1, L1f, L1D, L1Q, Lσf, Lσd, Lσq, LfD) -> ParkInductances

Calcula las inductancias de Park desde los parámetros del modelo de entrehierro.

Relaciones (ec. 1.7.7, 1.7.15, 1.7.16, 1.7.20):
    LG1 = (3/2)·L1,  LG2 = (3/2)·L2
    L1d = Lσ1 + LG1 + LG2
    L1q = Lσ1 + LG1 − LG2
    Lf1 = (3/2)·L1f,  LD1 = (3/2)·L1D,  LQ1 = (3/2)·L1Q
"""
function park_inductances(L1::Float64, L2::Float64, Lσ1::Float64,
                           L1f::Float64, L1D::Float64, L1Q::Float64,
                           Lσf::Float64, Lff_ag::Float64,
                           Lσd::Float64, LDD_ag::Float64,
                           Lσq::Float64, LQQ_ag::Float64,
                           LfD::Float64)
    LG1 = (3/2) * L1
    LG2 = (3/2) * L2
    return ParkInductances(
        Lσ1 + LG1 + LG2,          # L1d
        Lσ1 + LG1 - LG2,          # L1q
        Lσf + Lff_ag,              # Lf  (inductancia propia campo)
        Lσd + LDD_ag,              # LD
        Lσq + LQQ_ag,              # LQ
        (3/2) * L1f,               # Lf1
        (3/2) * L1D,               # LD1
        (3/2) * L1Q,               # LQ1
        LfD                        # LfD
    )
end

# -----------------------------------------------------------------------------
# Flujos de enlace en coordenadas dq (ec. 1.7.11 y análogas para rotor)
# ψ1d = L1d·i1d + Lf1·if + LD1·iD
# ψ1q = L1q·i1q + LQ1·iQ
# ψf  = Lf·if   + Lf1·i1d + LfD·iD
# ψD  = LD·iD   + LD1·i1d + LfD·if
# ψQ  = LQ·iQ   + LQ1·i1q
# -----------------------------------------------------------------------------

"""
    park_flux_linkages(i1d, i1q, if_, iD, iQ, li) -> (ψ1d, ψ1q, ψf, ψD, ψQ)

Calcula los enlaces de flujo en coordenadas dq para los cinco devanados
(estator d, estator q, campo, amortiguador D, amortiguador Q).
"""
function park_flux_linkages(i1d::Float64, i1q::Float64,
                             if_::Float64, iD::Float64, iQ::Float64,
                             li::ParkInductances)
    ψ1d = li.L1d * i1d + li.Lf1 * if_ + li.LD1 * iD
    ψ1q = li.L1q * i1q + li.LQ1 * iQ
    ψf  = li.Lf  * if_ + li.Lf1 * i1d + li.LfD * iD
    ψD  = li.LD  * iD  + li.LD1 * i1d + li.LfD * if_
    ψQ  = li.LQ  * iQ  + li.LQ1 * i1q
    return ψ1d, ψ1q, ψf, ψD, ψQ
end

# -----------------------------------------------------------------------------
# Ecuaciones de Park en forma de estado: dψ/dt
# Reescribiendo (ec. 1.7.13–1.7.14, 1.7.17–1.7.19) en forma de flujos:
#   dψ1d/dt = v1d − R1·i1d + ω·ψ1q
#   dψ1q/dt = v1q − R1·i1q − ω·ψ1d
#   dψf/dt  = vf  − Rf·if
#   dψD/dt  = vD  − RD·iD
#   dψQ/dt  = vQ  − RQ·iQ
# -----------------------------------------------------------------------------

"""
    park_flux_derivatives(v1d, v1q, vf, vD, vQ,
                          i1d, i1q, if_, iD, iQ,
                          ψ1d, ψ1q, R1, Rf, RD, RQ, ω)
    -> (dψ1d, dψ1q, dψf, dψD, dψQ)

Derivadas temporales de los enlaces de flujo en coordenadas dq.
ω = dγ/dt es la velocidad eléctrica del rotor [rad_elec/s].

Las tensiones de velocidad ±ω·ψ en los ejes d y q tienen origen en el
término j·ω·ψs,r de la ecuación compleja (ec. 1.7.10).
"""
function park_flux_derivatives(v1d::Float64, v1q::Float64,
                                vf::Float64,  vD::Float64,  vQ::Float64,
                                i1d::Float64, i1q::Float64,
                                if_::Float64, iD::Float64,  iQ::Float64,
                                ψ1d::Float64, ψ1q::Float64,
                                R1::Float64,  Rf::Float64,
                                RD::Float64,  RQ::Float64,
                                ω::Float64)
    dψ1d = v1d - R1 * i1d + ω * ψ1q
    dψ1q = v1q - R1 * i1q - ω * ψ1d
    dψf  = vf  - Rf * if_
    dψD  = vD  - RD * iD
    dψQ  = vQ  - RQ * iQ
    return dψ1d, dψ1q, dψf, dψD, dψQ
end

# -----------------------------------------------------------------------------
# Inversión del sistema de flujos para obtener corrientes
# Sistema lineal 5×5 desacoplado en dos subsistemas:
#   Eje d (3×3): [ψ1d, ψf, ψD] = Md · [i1d, if, iD]
#   Eje q (2×2): [ψ1q, ψQ]     = Mq · [i1q, iQ]
# -----------------------------------------------------------------------------

"""
    park_currents_from_flux(ψ1d, ψ1q, ψf, ψD, ψQ, li) -> (i1d, i1q, if_, iD, iQ)

Obtiene corrientes dq invirtiendo el sistema de inductancias.
Los subsistemas d y q están desacoplados.
"""
function park_currents_from_flux(ψ1d::Float64, ψ1q::Float64,
                                  ψf::Float64,  ψD::Float64, ψQ::Float64,
                                  li::ParkInductances)
    # Submatriz eje d: Md = [L1d Lf1 LD1; Lf1 Lf LfD; LD1 LfD LD]
    Md = @SMatrix [li.L1d  li.Lf1  li.LD1;
                   li.Lf1  li.Lf   li.LfD;
                   li.LD1  li.LfD  li.LD ]
    id_vec = Md \ SVector(ψ1d, ψf, ψD)

    # Submatriz eje q: Mq = [L1q LQ1; LQ1 LQ]
    Mq = @SMatrix [li.L1q  li.LQ1;
                   li.LQ1  li.LQ ]
    iq_vec = Mq \ SVector(ψ1q, ψQ)

    return id_vec[1], iq_vec[1], id_vec[2], id_vec[3], iq_vec[2]
end
