include("loadSPP.jl")
include("getfname.jl")

using LinearAlgebra

cd("$(homedir())/github/metaheuristiques")
#test ouvrir un fichier d'instance
fname = "DM1/Data/pb_100rnd0100.dat"
cost, matrix = loadSPP(fname)

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
    cost_restants = copy(cost)
    objective::Int = 0
    while length(sous_ensembles_restants) > 0
        sous_ensemble_candidat= zeros(Int,dim_matrix[2])
        max::Float64 = 0.0
        indice_choisi::Int = 0
        evaluation = Float64[]
        sizehint!(evaluation, dim_matrix[2])
#=
pour chaque sous ensemble restant, on prend celui qui maximise
cost/card(sous_ensemble) et qui est "disjoint" de la solution qu'on est en train de construire
=#
        for compteur in 1:length(sous_ensembles_restants)
            evaluation = push!(evaluation, cost_restants[compteur] / sum(matrix[sous_ensembles_restants[compteur],:]))
        end

        indice_choisi = argmax(evaluation)
        if verif_contrainte(x_0,matrix[sous_ensembles_restants[indice_choisi],:])
        x_0[indice_choisi] = 1
        end

        sous_ensembles_restants = deleteat!(sous_ensembles_restants,indice_choisi)
        cost_restants = deleteat!(cost_restants, indice_choisi)

        for compteur in sous_ensembles_restants
            if dot(x_0,matrix[compteur,:])>=1
                remove!(sous_ensembles_restants,compteur)
            end
        end
        println(sous_ensembles_restants)
    end
    objective = dot(x_0, cost)
    return x_0,objective
end

println(greedy_construction(cost,matrix))
