include("loadSPP.jl")
include("getfname.jl")

using LinearAlgebra

cd("$(homedir())/github/metaheuristiques")
#test ouvrir un fichier d'instance
fname = "DM1/Data/pb_100rnd0100.dat"
cost, matrix = loadSPP(fname)
# println(cost)
#println(matrix[1,:])

function verif_contrainte(vec1,vec2)
    return dot(vec1,vec2) < 1
end

#construction gloutonne d'une solution x_0 admissible et "de bonne facture"
#ne fonctionne pas
function greedy_construction(cost,matrix)
    dim_matrix = size(matrix)
    x_0 = zeros(Int, dim_matrix[2])

    cost_sorted = unique(sort(cost,rev=true)) #on trie les poids par ordre décroissant et on supprime les doublons

    #pour c
    for i in cost_sorted
        for j in 1:dim_matrix[1]
            if matrix[j,i] == 1 && verif_contrainte(x_0,matrix[j,:])
                x_0 = x_0 + matrix[j,:]
            end
        end
    end
    return x_0
end

dot(cost,greedy_construction(cost,matrix))
