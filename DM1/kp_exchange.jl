include("loadSPP.jl")
include("getfname.jl")

using LinearAlgebra

function exchange(k, p, x)
    if x == 0  && k < 2
        x = 1
        k = k + 1
    elseif x == 1 && p < 1
        x = 0
        p = p + 1
    end
    return k, p, x
end

function inv_exchange(k, p, x)
    if x == 0
        x = 1
        p = p - 1
    else
        x = 0
        k = k - 1
    end
    return k, p, x
end

function kp_exchange21(cost, matrix, x)
    k = 0
    p = 0
    max = 0
end
