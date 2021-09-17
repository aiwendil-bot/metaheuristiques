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

function remove!(a, item)
    deleteat!(a, findall(x->x==item, a))
end


#construction gloutonne d'une solution x_0 admissible et "de bonne facture"
function greedy_construction(cost,matrix)
    dim_matrix = size(matrix)
    x_0 = zeros(Int, dim_matrix[2])
    sous_ensembles_restants = 1:dim_matrix[2]
    sous_ensembles_restants = convert(Vector,sous_ensembles_restants)
    while length(sous_ensembles_restants) > 0
        sous_ensemble_candidat= zeros(Int,dim_matrix[2])
        max::Float16 = 0.0
        indice_choisi::Int = 0
        for compteur in sous_ensembles_restants
            if verif_contrainte(sous_ensemble_candidat,x_0) && cost[compteur]/sum(matrix[compteur,:])>max
                sous_ensemble_candidat = matrix[compteur,:]
                max = cost[compteur]/sum(matrix[compteur,:])
                indice_choisi = compteur
            end
        end
        sous_ensembles_restants = remove!(sous_ensembles_restants,indice_choisi)
        x_0 = x_0 + sous_ensemble_candidat
        for compteur in sous_ensembles_restants
            if dot(x_0,matrix[compteur,:])==1
                remove!(sous_ensembles_restants,compteur)
            end
        end
    end
    return x_0
end
println(greedy_construction(cost,matrix))
dot(cost,greedy_construction(cost,matrix))
