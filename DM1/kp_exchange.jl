include("loadSPP.jl")
include("getfname.jl")
include("greedyintelligent.jl")

using LinearAlgebra

fname = "DM1/Data/pb_100rnd0300.dat"
cost, matrix = loadSPP(fname)
x = greedy_intelligent(cost, matrix)
println(dot(cost, greedy_intelligent(cost, matrix)))

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

function kp_exchange21(cost, matrix, x0)
    k = 0
    p = 0
    x = x0
    max = dot(cost, x)
    better_x = Vector{Int64}(undef, length(x))
    for i in 1:length(x)
        k, p, x[i] = exchange(k, p, x[i])
        # println(x) # Test
        # println(dot(matrix[i,:], x)) # Test
        if dot(matrix[i,:], x) <= 1
            for j in i + i + 1:length(x)
                k, p, x[j] = exchange(k, p, x[j])
                # println(x) # Test
                # println(dot(matrix[j,:], x)) # Test
                if dot(matrix[j,:], x) <= 1
                    for h in j + 1:length(x)
                        k, p, x[h] = exchange(k, p, x[h])
                        eval = dot(cost, x)
                        if eval > max && dot(matrix[h,:], x) <= 1
                            max = eval
                            better_x = x
                        end
                        k, p, x[h] = inv_exchange(k, p, x[h])
                    end
                end
                k, p, x[j] = inv_exchange(k, p, x[j])
            end
        end
        k, p, x[i] = inv_exchange(k, p, x[i])
    end
    return better_x, max
end

kp_exchange21(cost, matrix, x)
