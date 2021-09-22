include("loadSPP.jl")
include("getfname.jl")

using LinearAlgebra

fname = "DM1/Data/pb_100rnd0100.dat"
cost, matrix = loadSPP(fname)

function greedy_construction(cost, matrix)
    copie_matrix = copy(matrix)
    x0 = zeros(size(matrix)[2])
    subset_restant = ones(size(matrix)[1])
    ratio = Vector{Float64}(undef, size(matrix)[2])
    println(size(ratio))
    while sum(subset_restant) > 0
        somme = sum(copie_matrix, dims=1)
        for i in 1:size(copie_matrix)[2]
            if somme[i] != 0
                ratio[i] = cost[i] / somme[i]
            end
        end
        indice_ratio_max = argmax(ratio)
        x0[indice_ratio_max] = 1
        subset_restant[indice_ratio_max] = 0
        for i in copie_matrix[indice_ratio_max, :]
            i = 0
        end
        for i in 1:size(matrix)[2]
            if dot(copie_matrix[i,:], x0) >= 1
                subset_restant[i] = 0
                for j in copie_matrix[i, :]
                    j = 0
                end
            end
        end
    end
    return x0, dot(x0, cost)
end

greedy_construction(cost, matrix)
