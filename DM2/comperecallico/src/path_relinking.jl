include("codeDM2.jl")
include("grasp.jl")

function grasppr(C,A,α,nb_iter,max_elite)
	ensemble_z_max = Vector{Int64}(undef,nb_iter)
	pool_elite = Dict()
	sizehint!(pool_elite,max_elite)
	x_max = Vector{Int64}(undef,length(C))
	z_max::Int64 = 0
	@distributed for i in 1:nb_iter
		x,z_greedy = greedy_randomized_construction(C,A,α)
		x,z = simple_descent(C,A,x,z_greedy)
		if i >= 2
			keys_elite = collect(keys(pool_elite))
			keys_elite = shuffle(keys_elite)
			y = keys_elite[1:rand(1:length(keys_elite))]
			for elite in y
				x_t = dot(C,x) > dot(C,elite) ? x : elite
				x_s = dot(C,x) < dot(C,elite) ? x : elite
				x_p, z_p = path_relinking(C, A, x_s, x_t)
				if length(pool_elite) < max_elite
					pool_elite[x_p] = z_p
				else
					minkey,minval = find_min_key(pool_elite)
					if z_p > minval
						delete!(pool_elite,minkey)
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
	return find_max_key(pool_elite),ensemble_z_max
end

function path_relinking(cost,matrix, x_s, x_t)
	symdiff = xor.(x_s, x_t)
	findall_symdiff = findall(x->x==1, symdiff) #indices des éléments qui ne sont pas dans l'intersection
	z_max = max(dot(cost,x_s),dot(cost,x_t))
	z_init = z_max
	x_max = dot(x_s,cost) > dot(x_t,cost) ? x_s : x_t
	x = copy(x_s)
	z_s = dot(cost,x_s)
	i=1
	while length(findall_symdiff) > 0

		valeurs_flips = Dict{Int64, Int64}()
		sizehint!(valeurs_flips, length(findall_symdiff))
		# indice_flip = rand(1:length(findall_symdiff)) # Passer de random à une évaluation
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
		#if max_key != 0
		x[max_key] = x[max_key] == 0 ? 1 : 0
		#end
		symdiff = xor.(x, x_t)
		findall_symdiff = findall(x->x==1, symdiff)
		if isAdmissible(cost,matrix,x)
			if max_val > z_init #solution "prometteuse" => autre critère ? z > 0.9 * z_max etc
				x, max_val = deepest_descent(cost,matrix,x,max_val)
			end

			if max_val > z_max
				z_max = max_val
				x_max = x
			end
		end
	end
	return x_max,z_max
end

function path_relinking_random(cost,matrix, x_s, x_t, i_max)
	symdiff = xor.(x_s, x_t)
	findall_symdiff = findall(x->x==1, symdiff) #indices des éléments qui ne sont pas dans l'intersection
	z_max = max(dot(cost,x_s),dot(cost,x_t))
	z_ini = z_max
	x_max = dot(x_s,cost) > dot(x_t,cost) ? x_s : x_t
	x = copy(x_s)
	while length(findall_symdiff) > 0# ?
		indice_flip = rand(1:length(findall_symdiff))
		x[findall_symdiff[indice_flip]] = x[findall_symdiff[indice_flip]] == 0 ? 1 : 0
		if est_admissible(x,matrix,findall_symdiff[indice_flip])
			z = dot(cost,x)
			if z > z_max #solution "prometteuse" => autre critère ? z > 0.9 * z_max etc
				x, z = deepest_descent(cost,matrix,x,z)
			end
			symdiff = xor.(x, x_t)
			findall_symdiff = findall(x->x==1, symdiff)
			if z > z_max
				z_max = z
				x_max = copy(x)
			end
		else
			x[findall_symdiff[indice_flip]] = x[findall_symdiff[indice_flip]] == 0 ? 1 : 0
		end
	end
	return x_max,z_max
end
