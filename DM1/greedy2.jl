include("loadSPP.jl")
include("getfname.jl")

using LinearAlgebra

fname = "DM1/Data/pb_100rnd0100.dat"
cost, matrix = loadSPP(fname)

function greedy_construction(cost, matrix)
    copie_matrix = copy(matrix)
    x0 = zeros(size(matrix)[1])
    subset_restant = ones(size(matrix)[1])
    while sum(subset_restant) > 0
        somme = sum(copie_matrix, dims=2)
        ratio = cost / somme
        indice_ratio_max = argmax(ratio)
        x0[indice_ratio_max] = 1
        subset_restant[indice_ratio_max] = 0
        for i in copie_matrix[indice_ratio_max, :]
            i = 0
        end
    end
    return x0, dot(x0, cost)
end

# greedy_construction(cost, matrix)
