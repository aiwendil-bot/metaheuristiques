include("loadSPP.jl")
include("getfname.jl")

using LinearAlgebra

cd("$(homedir())/github/metaheuristiques")

fname = "DM1/Data/pb_100rnd0100.dat"
cost, matrix = loadSPP(fname)

#renvoie un vecteur de vecteurs, de taille nb_variables
#chaque vecteur contient les indices des contraintes auxquelles la variable est liée

function vect_variables(matrix)
    result = [] #moche mais fuck la déclaration de vecteurs de vecteurs dans Julia
    sizehint!(result, size(matrix,2) )
    for i in 1:size(matrix,2)
        j = 1
        vect = Vector{Int}(undef,sum(matrix[:,i]))
        for k in 1:size(matrix,1)
            if matrix[k,i]==1
                vect[j] = k
                j = j + 1
            end
        end
    push!(result, vect)
    end
    return result
end

#renvoie un vecteur de vecteurs, de taille nb_contraintes
#chaque vecteur contient les indices des variable pour lesquelles la contrainte est liée

function vect_contraintes(matrix)
    result = [] #moche mais fuck la déclaration de vecteurs de vecteurs dans Julia
    sizehint!(result, size(matrix,1) )
    for i in 1:size(matrix,1)
        j = 1
        vect = Vector{Int}(undef,sum(matrix[i,:]))
        for k in 1:size(matrix,2)
            if matrix[i,k]==1
                vect[j] = k
                j = j + 1
            end
        end
    push!(result, vect)
    end
    return result
end
#calcul les évaluations retourne l'indice de la meilleure
function indice_meilleure_evaluation(cost,vect_variables,subset_restant,cost_restant)
	evaluations = zeros(Float64,length(cost))
	for i in 1:length(cost)
		somme = 0
		for j in 1:size(vect_variables[i],1)
			somme = somme + cost_restant[vect_variables[i][j]]
		end
		if somme != 0
			evaluations[i] = cost[i] / somme
		end
	end
	max,indice_max = 0, 0
	for k in 1:size(subset_restant,1)
		if evaluations[k]*subset_restant[k] > max
			indice_max = k
			max = evaluations[k]
		end
	end
	return indice_max
end

function greedy_intelligent(cost, matrix)
	vect_liaison_contraintes = vect_contraintes(matrix)
	vect_liaison_variables = vect_variables(matrix)

	x_0 = zeros(Int,size(vect_liaison_variables,1)) #on initie la solution avec des 0

	sous_ensembles_restants = ones(Int,size(vect_liaison_variables,1)) # si 0 alors sous ensemble vu
	cost_restants = ones(Int,size(vect_liaison_contraintes,1)) # si 0 alors coût utilisé

	while sum(sous_ensembles_restants) > 0 #tant qu'il reste des sous ensembles non vus
		#on calcule les nouvelles évaluations :
		indice_max = indice_meilleure_evaluation(cost,vect_liaison_variables,sous_ensembles_restants,cost_restants)
		#le cas indice_max==0 peut survenir s'il ne reste que des sous_ensembles avec des coûts égaux à 0, auquel cas on les prend
		if indice_max == 0
			n = find(x->x==1,sous_ensembles_restants)
			for k in n
				x_0[k] = 1
				sous_ensembles_restants[k] = 0
			end
		else
			x_0[indice_max] = 1
			sous_ensembles_restants[indice_max] = 0
			for i in vect_liaison_variables[indice_max]
				if cost_restants[i] != 0
					for k in vect_liaison_contraintes[i]
						if k != indice_max && sous_ensembles_restants[k] == 1
							x_0[k] = 0
							sous_ensembles_restants[k] = 0
						end
					end
					cost_restants[i] = 0
				end
			end
		end
	end
	return x_0
end

println(greedy_intelligent(cost, matrix))
dot(cost,greedy_intelligent(cost,matrix))
