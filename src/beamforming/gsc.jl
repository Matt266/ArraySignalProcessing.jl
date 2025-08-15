"""
References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002.
"""
struct GeneralizedSidelobeCanceler{UT, WQT, WAT, BT} <: AbstractBeamformer where {UT<:AbstractUpdateMethod}
    UpdateMethod:: UT
    WQ::WQT
    WA::WAT
    B::BT
end

"""
Non-adaptive filtering of the input data matrix X
"""
function process(gsc::GeneralizedSidelobeCanceler, X)
    # GSC outputs
    Yc = gsc.WQ'*X
    Z = gsc.B'*X
    Yb = gsc.WA'*Z
    Y = Yc - Yb
    return Y
end

"""
Adaptive filtering of the input data matrix X. Returns X, mutates the GSC

References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002.
"""
function process!(gsc::GeneralizedSidelobeCanceler{<:RLS}, X)
    # short name so formulas can be read better
    μ = gsc.UpdateMethod.μ
    P = gsc.UpdateMethod.P

    # GSC outputs
    Yc = gsc.WQ'*X
    Z = gsc.B'*X
    Yb = gsc.WA'*Z
    Y = Yc - Yb

    # RLS update
    #gz = 1/μ*P*Z/(1+1/μ*Z'*P*Z) #(7.199)
    G = P*Z/(μ*I+Z'*P*Z)

    # P  1/μ*P - 1/μ*gz*Z'*P #(7.200)
    gsc.UpdateMethod.P .= (P-G*Z'*P)/μ

    # Beamformer weight update
    # wa += gz*conj(Y) #(7.204)
    gsc.WA .+= G*Y'

    return Y
end

"""
References:
-----------
H. L. Van Trees, Optimum array processing. Nashville, TN: John Wiley & Sons, 2002.
"""
function process!(gsc::GeneralizedSidelobeCanceler{<:LMS}, X)
    # short name so formulas can be read better
    μ = gsc.UpdateMethod.μ

    # GSC outputs
    Yc = gsc.WQ'*X
    Z = gsc.B'*X
    Yb = gsc.WA'*Z
    Y = Yc - Yb

    # Beamformer weight update
    # (7.425), (7.426)
    gsc.WA .+= μ*Z*Y'

    return Y
end