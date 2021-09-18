include("loadSPP.jl")
include("getfname.jl")

using LinearAlgebra

cd("$(homedir())/github/metaheuristiques")

fname = "DM1/Data/pb_100rnd0100.dat"
cost, matrix = loadSPP(fname)

function verif_contrainte(vec1,vec2)
    return dot(vec1,vec2) < 1
end

function remove!(tab, item)
    deleteat!(tab, findall(x->x==item, tab))
end


#construction gloutonne d'une solution x_0 admissible et "de bonne facture"
function greedy_construction(cost,matrix)
    dim_matrix = size(matrix)
    x_0 = zeros(Int, dim_matrix[2])
    sous_ensembles_restants = 1:dim_matrix[2]
    sous_ensembles_restants = convert(Vector,sous_ensembles_restants)
    objective::Int = 0
    while length(sous_ensembles_restants) > 0
        evaluations = Float64[]
        sizehint!(evaluations, dim_matrix[2])
#=
pour chaque sous ensemble restant, on prend celui qui maximise
cost/card(sous_ensemble) et qui est "disjoint" de la solution qu'on est en train de construire
=#
        for compteur in 1:length(sous_ensembles_restants)
            evaluations = push!(evaluations, cost[sous_ensembles_restants[compteur]]/ sum(matrix[sous_ensembles_restants[compteur],:]))
        end

        indice_choisi = argmax(evaluations)

        x_0[sous_ensembles_restants[indice_choisi]] = 1

#le sous ensemble choisi et son coût associé sont supprimés
        sous_ensembles_restants = deleteat!(sous_ensembles_restants,indice_choisi)

#on identifie et supprime les sous_ensembles qui ne pourront plus être pris

        for compteur in sous_ensembles_restants
            if !verif_contrainte(x_0,matrix[compteur,:])
                sous_ensembles_restants = remove!(sous_ensembles_restants, compteur)
            end
        end
    end
    objective = dot(x_0, cost)
    return x_0,objective
end
println(sum(greedy_construction(cost,matrix)[1]))
println(greedy_construction(cost,matrix))
