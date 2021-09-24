include("loadSPP.jl")
include("getfname.jl")

using LinearAlgebra

cd("$(homedir())/github/metaheuristiques")

fname = "DM1/Data/pb_100rnd0100.dat"
cost, matrix = loadSPP(fname)

# Pour savoir si c'est disjoint
function verif_contrainte(vec1, vec2)
    return dot(vec1, vec2) < 1
end

# Raccourci
function remove!(tab, item)
    deleteat!(tab, findall(x -> x == item, tab))
end


# construction gloutonne d'une solution x_0 admissible et "de bonne facture"
function greedy_construction(cost, matrix)
    somme = sum(matrix, dims=2)
    dim_matrix = size(matrix, 2)
    x_0 = zeros(Int, dim_matrix)
    sous_ensembles_restants = [x for x in 1:dim_matrix]
    println(sous_ensembles_restants)
    # sous_ensembles_restants = convert(Vector, sous_ensembles_restants)
    objective::Int = 0
    while length(sous_ensembles_restants) > 0
        evaluations = Float64[]
        sizehint!(evaluations, dim_matrix)
#=
pour chaque sous ensemble restant, on prend celui qui maximise
cost/card(sous_ensemble) =#
        for compteur in 1:length(sous_ensembles_restants)
            push!(evaluations, cost[sous_ensembles_restants[compteur]] / somme[sous_ensembles_restants[compteur]])
        end
        println(evaluations)

        indice_choisi = argmax(evaluations) # Amélioration possible : si 2 argmax identique, vérifier la ligne avec le moins de 1

        x_0[sous_ensembles_restants[indice_choisi]] = 1

# le sous ensemble choisi est supprimé
        deleteat!(sous_ensembles_restants, indice_choisi)

# on identifie et supprime les sous_ensembles qui ne pourront plus être pris

        for compteur in sous_ensembles_restants
            if !verif_contrainte(x_0, matrix[compteur,:])
                remove!(sous_ensembles_restants, compteur)
            end
        end
        # println(sous_ensembles_restants)
    end
    objective = dot(x_0, cost)
    return x_0, objective
end
# println(sum(greedy_construction(cost, matrix)[1]))
println(greedy_construction(cost, matrix))
