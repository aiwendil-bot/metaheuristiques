include("loadSPP.jl")
include("getfname.jl")
include("greedyintelligent.jl")

using LinearAlgebra

fname = "DM1/Data/pb_200rnd0900.dat"
cost, matrix = loadSPP(fname)
x = greedy_intelligent(cost, matrix)
z = dot(cost, greedy_intelligent(cost, matrix))
println(z)
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

#kp_exchange21(cost, matrix, x)


function kp_exchange01(cost, matrix, x0, z)
    x = copy(x0)
    z_ameliore = copy(z)
    max::Int64 = copy(z)
    for i in 1:length(x)
        if x0[i] == 0
            x[i] = 1
            if est_admissible(x,matrix)
                z_ameliore = z + cost[i]
                if z_ameliore > max
                    max = z_ameliore
                    x0 = copy(x)
                    println("x")
                end
            end
            x[i] =0
        end
    end
    return max
end

function est_admissible(x,matrix)
    for compteur in 1:size(matrix,1)
        if dot(matrix[compteur,:],x) >1
            return false
        end
    end
    return true
end

#kp_exchange01(cost, matrix, x,z)

function kp_exchange11(cost, matrix, x0, z)
    x = copy(x0)
    z_ameliore = copy(z)
    max::Int64 = copy(z)
    for i in 1:length(x)
        if x0[i] == 0
            for j in 1:length(x)
                if x0[j] == 1
                    x[i],x[j] = 1,0
                    if est_admissible(x,matrix)
                        z_ameliore = z + cost[i] - cost[j]
                        if z_ameliore > max
                            max = z_ameliore
                            x0 = copy(x)
                            println("xx")
                        end
                    end
                    x0[i],x[j] = 0,1
                end
            end
        end
    end
    return max
end
@time kp_exchange11(cost, matrix, x,z)
