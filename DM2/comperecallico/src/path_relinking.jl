include("codeDM2.jl")
include("grasp.jl")

function grasppr(C, liaisons_contraintes, liaisons_variables, α, nb_iter, max_elite)
	ensemble_z_max = Vector{Int64}(undef, nb_iter)
	pool_elite = Dict()
	sizehint!(pool_elite, max_elite)
	x_max = Vector{Int64}(undef, length(C))
	z_max::Int64 = 0
	t1::Float64,t2::Float64,t::Float64,check::Bool = 0.0,0.0,0.0,false

	for i in 1:nb_iter


		x, z_greedy = greedy_randomized_construction(C, liaisons_contraintes, liaisons_variables, α)
		x, z = simple_descent(C, liaisons_contraintes, liaisons_variables, x, z_greedy)
		if i >= 2
			t1 = time()
			keys_elite = collect(keys(pool_elite))
			keys_elite = shuffle(keys_elite)
			y = keys_elite[1:rand(1:length(keys_elite))]
			for elite in y
				x_t = dot(C, x) > dot(C, elite) ? x : elite
				x_s = dot(C, x) < dot(C, elite) ? x : elite
				x_p, z_p = path_relinking(C, liaisons_contraintes, liaisons_variables, x_s, x_t)

				if length(pool_elite) < max_elite
					pool_elite[x_p] = z_p
				else
					minkey, minval = find_min_key(pool_elite)
					if z_p > minval
						delete!(pool_elite, minkey)
						pool_elite[x_p] = z_p
					end
				end
				if z_p > z_max
					z_max = z_p
					x_max = x_p
				end
			end
		else
			pool_elite[x] = z
		end

		ensemble_z_max[i] = find_max_key(pool_elite)[2]
	end
	return find_max_key(pool_elite), dot(C,find_max_key(pool_elite)[1])
end

function path_relinking(cost, liaisons_contraintes, liaisons_variables, x_s, x_t)
	symdiff = xor.(x_s, x_t)
	findall_symdiff = findall(x -> x == 1, symdiff) # indices des éléments qui ne sont pas dans l'intersection
	z_max = max(dot(cost, x_s), dot(cost, x_t))
	z_init = z_max
	x_max = dot(x_s, cost) > dot(x_t, cost) ? x_s : x_t
	x = copy(x_s)
	z_s = dot(cost, x_s)
	i = 1
	while length(findall_symdiff) > 0
		valeurs_flips = Dict{Int64,Int64}()
		sizehint!(valeurs_flips, length(findall_symdiff))
		for i in 1:length(findall_symdiff)
			x[findall_symdiff[i]] = x[findall_symdiff[i]] == 0 ? 1 : 0
			z = x[findall_symdiff[i]] == 1 ? z_s + cost[findall_symdiff[i]] : z_s - cost[findall_symdiff[i]]
			valeurs_flips[findall_symdiff[i]] = z
			x[findall_symdiff[i]] = x[findall_symdiff[i]] == 0 ? 1 : 0
        end
		max_key, max_val = 0, 0
		for i in keys(valeurs_flips)
			if valeurs_flips[i] > max_val
				max_key = i
				max_val = valeurs_flips[i]
			end
		end
		if max_key != 0
			x[max_key] = x[max_key] == 0 ? 1 : 0
			symdiff = xor.(x, x_t)

		end
		if max_val == 0
			symdiff = zeros(length(x))
		end
		findall_symdiff = findall(x -> x == 1, symdiff)
		admis::Bool = true
		j = 1
		while j <= length(cost) && admis == true
			admis = est_admissible(cost, liaisons_contraintes, liaisons_variables, j)
			j += 1
		end
		if admis == true
			if max_val > z_init # solution "prometteuse" => autre critère ? z > 0.9 * z_max etc
				x, max_val = deepest_descent(cost, liaisons_contraintes, liaisons_variables, x, max_val)
			end
			if max_val > z_max
				z_max = max_val
				x_max = x
			end
		end
	end
	return x_max, z_max
end

function path_relinking_random(cost, liaisons_contraintes, liaisons_variables, x_s, x_t)
	symdiff = xor.(x_s, x_t)
	findall_symdiff = findall(x -> x == 1, symdiff) # indices des éléments qui ne sont pas dans l'intersection
	z_max = max(dot(cost, x_s), dot(cost, x_t))
	z_ini = z_max
	x_max = dot(x_s, cost) > dot(x_t, cost) ? x_s : x_t
	x = copy(x_s)
	while length(findall_symdiff) > 0
		indice_flip = rand(1:length(findall_symdiff))
		x[findall_symdiff[indice_flip]] = x[findall_symdiff[indice_flip]] == 0 ? 1 : 0

			z = dot(cost, x)
			if z > z_max # solution "prometteuse" => autre critère ? z > 0.9 * z_max etc
				if est_admissible(x, liaisons_contraintes, liaisons_variables, findall_symdiff[indice_flip])
					x, z = deepest_descent(cost, liaisons_contraintes, liaisons_variables, x, z)
				end
			end
			symdiff = xor.(x, x_t)
			findall_symdiff = findall(x -> x == 1, symdiff)
			if z > z_max
				z_max = z
				x_max = copy(x)
			end
	end
	return x_max, z_max
end
