# --------------------------------------------------------------------------- #
# Perform a numerical experiment (with a fake version of GRASP-SPP)
#= @everywhere include("grasp.jl")
@everywhere include("path_relinking.jl")
@everywhere include("reactive_grasp.jl")
@everywhere include("../../../libSPP/librarySPP.jl") =#
using PyPlot
pygui(true)

#GRASP AVEC OU SANS REACTIVE GRASP
function graspSPP(cost,liaisons_contraintes,liaisons_variables, alpha, nbIterationGrasp)

    zconstruction = zeros(Int64, nbIterationGrasp)
    zamelioration = zeros(Int64, nbIterationGrasp)
    zbest = zeros(Int64, nbIterationGrasp)
    zbetter = 0

    for i = 1:nbIterationGrasp
        zconstruction[i], zamelioration[i] = grasp_v2(cost, liaisons_contraintes, liaisons_variables, alpha)[2], grasp_v2(cost, liaisons_contraintes, liaisons_variables, alpha)[4]# # livrable du DM2
        zbetter = max(zbetter, zamelioration[i])
        zbest[i] = zbetter
    end
    return zconstruction, zamelioration, zbest
end

function graspSPP_Reactive(cost,liaisons_contraintes,liaisons_variables, vector_α, nbIterationGrasp,N_α)

    zconstruction = zeros(Int64, nbIterationGrasp)
    zamelioration = zeros(Int64, nbIterationGrasp)
    zbest = zeros(Int64, nbIterationGrasp)
    zbetter = 0

    nb_α = length(vector_α)
    probas_α = 1 / nb_α .* ones(Float64, nb_α)
    probas_α_cumulees = zeros(Float64, nb_α)
    q = zeros(Float64, nb_α) # notation cours
    compteur_chaque_α = zeros(Int64, nb_α)
    z_best::Int64, z_worst::Int64, z_average = 0, typemax(Int64), Vector{Float64}(undef, nb_α)
    x_max = zeros(Int64, length(cost))
    z_glouton::Int64 = 0

    for i in 1:nbIterationGrasp
        s_prob = 0.0
        for j in 1:nb_α
            s_prob += probas_α[j]
            probas_α_cumulees[j] = s_prob
        end
        probas_α_cumulees[nb_α] = 1.0 # arrondis
            a = rand()
            indice_α = findfirst(x -> x >= a, probas_α_cumulees)
            x, z_glouton = greedy_randomized_construction(cost, liaisons_contraintes, liaisons_variables, vector_α[indice_α])
            x, z = simple_descent(cost, liaisons_contraintes, liaisons_variables, x, z_glouton)
            if z > z_best
                z_best = z
                x_max = x
            end
            if z < z_worst
                z_worst = z
            end
            z_average[indice_α] = (compteur_chaque_α[indice_α] * z_average[indice_α] + z) / (compteur_chaque_α[indice_α] + 1)
            compteur_chaque_α[indice_α] += 1

        if i % N_α == 0
            for k in 1:nb_α
                q[k] = (z_average[k] - z_worst) / (z_best - z_worst)
            end
            probas_α = q / sum(q)
        end
        zconstruction[i], zamelioration[i] = z_glouton, z
        zbetter = max(zbetter, zamelioration[i])
        zbest[i] = zbetter
    end
    return zconstruction, zamelioration, zbest
end

#GRASP AVEC PATH RELINKING
function graspSPP_PR(cost,liaisons_contraintes,liaisons_variables, alpha, nbIterationGrasp,max_elite)

    zconstruction = zeros(Int64, nbIterationGrasp)
    zamelioration = zeros(Int64, nbIterationGrasp)
    zbest = zeros(Int64, nbIterationGrasp)
    zbetter = 0

    ensemble_z_max = Vector{Int64}(undef, nbIterationGrasp)
    pool_elite = Dict()
    sizehint!(pool_elite, max_elite)
    x_max = Vector{Int64}(undef, length(cost))
    z_max::Int64 = 0
    z_p::Int64 = 0

    for i = 1:nbIterationGrasp
        x, z_greedy = greedy_randomized_construction(cost, liaisons_contraintes, liaisons_variables, alpha)
        x, z = simple_descent(cost, liaisons_contraintes, liaisons_variables, x, z_greedy)
        if i >= 2
            keys_elite = collect(keys(pool_elite))
            keys_elite = shuffle(keys_elite)
            y = keys_elite[1:rand(1:length(keys_elite))]
            for elite in y
                x_t = dot(cost, x) > dot(cost, elite) ? x : elite
                x_s = dot(cost, x) < dot(cost, elite) ? x : elite
                x_p, z_p = path_relinking(cost, liaisons_contraintes, liaisons_variables, x_s, x_t)
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
        zconstruction[i], zamelioration[i] = z_greedy, z_p
        zbetter = max(zbetter, zamelioration[i])
        zbest[i] = zbetter
    end
    #return find_max_key(pool_elite), ensemble_z_max

    return zconstruction, zamelioration, zbest
end


function plotRunGrasp(iname, zinit, zls, zbest)
    figure("Examen d'un run", figsize=(6, 6)) # Create a new figure
    title("GRASP-SPP | \$z_{Init}\$  \$z_{LS}\$  \$z_{Best}\$ | " * iname)
    xlabel("Itérations")
    ylabel("valeurs de z(x)")
    ylim(0, maximum(zbest) + 2)

    nPoint = length(zinit)
    x = collect(1:nPoint)
    xticks([1,convert(Int64, ceil(nPoint / 4)),convert(Int64, ceil(nPoint / 2)), convert(Int64, ceil(nPoint / 4 * 3)),nPoint])
    plot(x, zbest, linewidth=2.0, color="green", label="meilleures solutions")
    plot(x, zls, ls="", marker="^", ms=2, color="green", label="toutes solutions améliorées")
    plot(x, zinit, ls="", marker=".", ms=2, color="red", label="toutes solutions construites")
    vlines(x, zinit, zls, linewidth=0.5)
    legend(loc=4, fontsize="small")
end

function plotAnalyseGrasp(iname, x, zmoy, zmin, zmax)
    figure("bilan tous runs", figsize=(6, 6)) # Create a new figure
    title("GRASP-SPP | \$z_{min}\$  \$z_{moy}\$  \$z_{max}\$ | " * iname)
    xlabel("Itérations (pour nbRunGrasp)")
    ylabel("valeurs de z(x)")
    ylim(0, zmax[end] + 2)

    nPoint = length(x)
    intervalle = [reshape(zmoy, (1, nPoint)) - reshape(zmin, (1, nPoint)) ; reshape(zmax, (1, nPoint)) - reshape(zmoy, (1, nPoint))]
    xticks(x)
    errorbar(x, zmoy, intervalle, lw=1, color="black", label="zMin zMax")
    plot(x, zmoy, linestyle="-", marker="o", ms=4, color="green", label="zMoy")
    legend(loc=4, fontsize="small")
end

function plotCPUt(allfinstance, tmoy)
    figure("bilan CPUt tous runs", figsize=(6, 6)) # Create a new figure
    title("GRASP-SPP | tMoy")
    ylabel("CPUt moyen (s)")

    xticks(collect(1:length(allfinstance)), allfinstance, rotation=60, ha="right")
    margins(0.15)
    subplots_adjust(bottom=0.15, left=0.21)
    plot(collect(1:length(allfinstance)), tmoy, linestyle="--", lw=0.5, marker="o", ms=4, color="blue", label="tMoy")
    legend(loc=4, fontsize="small")
end


# Simulation d'une experimentation numérique  --------------------------

# Pkg.add("PyPlot") # Mandatory before the first use of this package
using PyPlot

function simulation()
    allfinstance      =  ["../../../Data/pb_100rnd0100.dat"]
    nbInstances       =  length(allfinstance)
    nbRunGrasp        =  30   # nombre de fois que la resolution GRASP est repetee
    nbIterationGrasp  =  200  # nombre d'iteration que compte une resolution GRASP
    nbDivisionRun     =  10   # nombre de division que compte une resolution GRASP

    zinit = zeros(Int64, nbIterationGrasp) # zero
    zls   = zeros(Int64, nbIterationGrasp) # zero
    zbest = zeros(Int64, nbIterationGrasp) # zero

    x     = zeros(Int64, nbDivisionRun)
    zmax  = Matrix{Int64}(undef, nbInstances, nbDivisionRun); zmax[:] .= typemin(Int64)  # -Inf entier
    zmoy  = zeros(Float64, nbInstances, nbDivisionRun) # zero
    zmin  = Matrix{Int64}(undef, nbInstances, nbDivisionRun) ; zmin[:] .= typemax(Int64)  # +Inf entier
    tmoy  = zeros(Float64, nbInstances)  # zero

    # calcule la valeur du pas pour les divisions
    for division = 1:nbDivisionRun
        x[division] = convert(Int64, ceil(nbIterationGrasp / nbDivisionRun * division))
    end

    println("Experimentation GRASP-SPP avec :")
    println("  nbInstances       = ", nbInstances)
    println("  nbRunGrasp        = ", nbRunGrasp)
    println("  nbIterationGrasp  = ", nbIterationGrasp)
    println("  nbDivisionRun     = ", nbDivisionRun)
    println(" ")
    cpt = 0
    cost, matrix = loadSPP(allfinstance[1])
    liaisons_contraintes = vect_contraintes(matrix)
    liaisons_variables = vect_variables(matrix)
    # run non comptabilise (afin de produire le code compile)
    zinit, zls, zbest = graspSPP(cost,liaisons_contraintes,liaisons_variables, 0.5, 1)

    for instance = 1:nbInstances
        # les instances sont traitees separement
        cost, matrix = loadSPP(allfinstance[instance])
        liaisons_contraintes = vect_contraintes(matrix)
        liaisons_variables = vect_variables(matrix)

        print("  ", allfinstance[instance], " : ")
        for runGrasp = 1:nbRunGrasp
            # une instance sera resolue nbrungrasp fois

            start = time() # demarre le compteur de temps
            alpha = 0.75

            zinit, zls, zbest = graspSPP(cost, liaisons_contraintes, liaisons_variables, alpha, nbIterationGrasp)
            tutilise = time() - start # arrete et releve le compteur de temps
            cpt += 1; print(cpt % 10)

            # mise a jour des resultats collectes
            for division = 1:nbDivisionRun
                zmax[instance,division] = max(zbest[x[division]], zmax[instance,division])
                zmin[instance,division] = min(zbest[x[division]], zmin[instance,division])
                zmoy[instance,division] =  zbest[x[division]] + zmoy[instance,division]
            end # division
            tmoy[instance] = tmoy[instance] + tutilise

        end # run
        for division = 1:nbDivisionRun
            zmoy[instance,division] =  zmoy[instance,division] /  nbRunGrasp
        end # division
        tmoy[instance] = tmoy[instance] / nbRunGrasp
        println(" ")

    end # instance

    # Pkg.add("PyPlot") # Mandatory before the first use of this package
    println(" ");println("  Graphiques de synthese")
#    using PyPlot
    instancenb = 1
    plotRunGrasp(allfinstance[instancenb], zinit, zls, zbest)
    plotAnalyseGrasp(allfinstance[instancenb], x, zmoy[instancenb,:], zmin[instancenb,:], zmax[instancenb,:])
    plotCPUt(allfinstance, tmoy)
end

function simulation_Reactive()
    allfinstance      =  ["../../../Data/pb_100rnd0100.dat"]
    nbInstances       =  length(allfinstance)
    nbRunGrasp        =  30   # nombre de fois que la resolution GRASP est repetee
    nbIterationGrasp  =  200  # nombre d'iteration que compte une resolution GRASP
    nbDivisionRun     =  10   # nombre de division que compte une resolution GRASP

    zinit = zeros(Int64, nbIterationGrasp) # zero
    zls   = zeros(Int64, nbIterationGrasp) # zero
    zbest = zeros(Int64, nbIterationGrasp) # zero

    x     = zeros(Int64, nbDivisionRun)
    zmax  = Matrix{Int64}(undef, nbInstances, nbDivisionRun); zmax[:] .= typemin(Int64)  # -Inf entier
    zmoy  = zeros(Float64, nbInstances, nbDivisionRun) # zero
    zmin  = Matrix{Int64}(undef, nbInstances, nbDivisionRun) ; zmin[:] .= typemax(Int64)  # +Inf entier
    tmoy  = zeros(Float64, nbInstances)  # zero

    # calcule la valeur du pas pour les divisions
    for division = 1:nbDivisionRun
        x[division] = convert(Int64, ceil(nbIterationGrasp / nbDivisionRun * division))
    end

    println("Experimentation Reactive-GRASP-SPP avec :")
    println("  nbInstances       = ", nbInstances)
    println("  nbRunGrasp        = ", nbRunGrasp)
    println("  nbIterationGrasp  = ", nbIterationGrasp)
    println("  nbDivisionRun     = ", nbDivisionRun)
    println(" ")
    cpt = 0
    cost, matrix = loadSPP(allfinstance[1])
    liaisons_contraintes = vect_contraintes(matrix)
    liaisons_variables = vect_variables(matrix)
    # run non comptabilise (afin de produire le code compile)
    zinit, zls, zbest = graspSPP(cost,liaisons_contraintes,liaisons_variables, 0.5, 1)

    for instance = 1:nbInstances
        # les instances sont traitees separement
        cost, matrix = loadSPP(allfinstance[instance])
        liaisons_contraintes = vect_contraintes(matrix)
        liaisons_variables = vect_variables(matrix)

        print("  ", allfinstance[instance], " : ")
        for runGrasp = 1:nbRunGrasp
            # une instance sera resolue nbrungrasp fois

            start = time() # demarre le compteur de temps
            alpha = 0.75
            vector_α = [0.2,0.5,0.7,0.9]
            N_α = 20

            zinit, zls, zbest = graspSPP_Reactive(cost,liaisons_contraintes,liaisons_variables, vector_α, nbIterationGrasp,N_α)
            tutilise = time() - start # arrete et releve le compteur de temps
            cpt += 1; print(cpt % 10)

            # mise a jour des resultats collectes
            for division = 1:nbDivisionRun
                zmax[instance,division] = max(zbest[x[division]], zmax[instance,division])
                zmin[instance,division] = min(zbest[x[division]], zmin[instance,division])
                zmoy[instance,division] =  zbest[x[division]] + zmoy[instance,division]
            end # division
            tmoy[instance] = tmoy[instance] + tutilise

        end # run
        for division = 1:nbDivisionRun
            zmoy[instance,division] =  zmoy[instance,division] /  nbRunGrasp
        end # division
        tmoy[instance] = tmoy[instance] / nbRunGrasp
        println(" ")

    end # instance

    # Pkg.add("PyPlot") # Mandatory before the first use of this package
    println(" ");println("  Graphiques de synthese")
#    using PyPlot
    instancenb = 1
    plotRunGrasp(allfinstance[instancenb], zinit, zls, zbest)
    plotAnalyseGrasp(allfinstance[instancenb], x, zmoy[instancenb,:], zmin[instancenb,:], zmax[instancenb,:])
    plotCPUt(allfinstance, tmoy)
end

function simulation_PR()
    allfinstance      =  ["../../../Data/pb_100rnd0100.dat"]
    nbInstances       =  length(allfinstance)
    nbRunGrasp        =  30   # nombre de fois que la resolution GRASP est repetee
    nbIterationGrasp  =  200  # nombre d'iteration que compte une resolution GRASP
    nbDivisionRun     =  10   # nombre de division que compte une resolution GRASP

    zinit = zeros(Int64, nbIterationGrasp) # zero
    zls   = zeros(Int64, nbIterationGrasp) # zero
    zbest = zeros(Int64, nbIterationGrasp) # zero

    x     = zeros(Int64, nbDivisionRun)
    zmax  = Matrix{Int64}(undef, nbInstances, nbDivisionRun); zmax[:] .= typemin(Int64)  # -Inf entier
    zmoy  = zeros(Float64, nbInstances, nbDivisionRun) # zero
    zmin  = Matrix{Int64}(undef, nbInstances, nbDivisionRun) ; zmin[:] .= typemax(Int64)  # +Inf entier
    tmoy  = zeros(Float64, nbInstances)  # zero

    # calcule la valeur du pas pour les divisions
    for division = 1:nbDivisionRun
        x[division] = convert(Int64, ceil(nbIterationGrasp / nbDivisionRun * division))
    end

    println("Experimentation GRASP-SPP et Path-relinking avec :")
    println("  nbInstances       = ", nbInstances)
    println("  nbRunGrasp        = ", nbRunGrasp)
    println("  nbIterationGrasp  = ", nbIterationGrasp)
    println("  nbDivisionRun     = ", nbDivisionRun)
    println(" ")
    cpt = 0
    cost, matrix = loadSPP(allfinstance[1])
    liaisons_contraintes = vect_contraintes(matrix)
    liaisons_variables = vect_variables(matrix)
    # run non comptabilise (afin de produire le code compile)
    zinit, zls, zbest = graspSPP_PR(cost,liaisons_contraintes,liaisons_variables, 0.5, 1,1)

    for instance = 1:nbInstances
        # les instances sont traitees separement
        cost, matrix = loadSPP(allfinstance[instance])
        liaisons_contraintes = vect_contraintes(matrix)
        liaisons_variables = vect_variables(matrix)

        print("  ", allfinstance[instance], " : ")
        for runGrasp = 1:nbRunGrasp
            # une instance sera resolue nbrungrasp fois

            start = time() # demarre le compteur de temps
            alpha = 0.75
            max_elite = 5

            zinit, zls, zbest = graspSPP_PR(cost, liaisons_contraintes, liaisons_variables, alpha, nbIterationGrasp,max_elite)
            tutilise = time() - start # arrete et releve le compteur de temps
            cpt += 1; print(cpt % 10)

            # mise a jour des resultats collectes
            for division = 1:nbDivisionRun
                zmax[instance,division] = max(zbest[x[division]], zmax[instance,division])
                zmin[instance,division] = min(zbest[x[division]], zmin[instance,division])
                zmoy[instance,division] =  zbest[x[division]] + zmoy[instance,division]
            end # division
            tmoy[instance] = tmoy[instance] + tutilise

        end # run
        for division = 1:nbDivisionRun
            zmoy[instance,division] =  zmoy[instance,division] /  nbRunGrasp
        end # division
        tmoy[instance] = tmoy[instance] / nbRunGrasp
        println(" ")

    end # instance

    # Pkg.add("PyPlot") # Mandatory before the first use of this package
    println(" ");println("  Graphiques de synthese")
#    using PyPlot
    instancenb = 1
    plotRunGrasp(allfinstance[instancenb], zinit, zls, zbest)
    plotAnalyseGrasp(allfinstance[instancenb], x, zmoy[instancenb,:], zmin[instancenb,:], zmax[instancenb,:])
    plotCPUt(allfinstance, tmoy)
end

function reactive_experiment(cost, liaisons_contraintes, liaisons_variables, vector_α, nb_iter, N_α)
    x = reactive_grasp(cost, liaisons_contraintes, liaisons_variables, vector_α, nb_iter, N_α)
    figure("test", figsize=(8, 8))
    title("Probabilité de chaque α pour $nb_iter itérations | refresh tous les $N_α")
    pie(x[3], labels=["α =" * string(i) for i in vector_α], normalize=true, autopct="%1.1f%%")
end

function comparaisons(C, liaisons_contraintes, liaisons_variables, α, nb_iter, max_elite)
    w = reactive_grasp(C, liaisons_contraintes, liaisons_variables, [0.2,0.4,0.6,0.8,0.9], nb_iter, 20)
    y = grasp(C, liaisons_contraintes, liaisons_variables, α, nb_iter)
    z = grasppr(C, liaisons_contraintes, liaisons_variables, α, nb_iter, max_elite)
    x = [i for i in 1:nb_iter]
    figure("comparaisons", figsize=(8, 8))
    title("Comparaison entre GRASP, Reactive-GRASP & Path-relinking | α = $α")
    xlabel("nombre d'itérations")
    ylabel("z_max")
    plot(x, w[4], linestyle="--", marker=".", color="green", label="Reactive-GRASP")
    plot(x, y[3], linestyle="--", marker=".", color="red", label="GRASP")
    plot(x, z[2], linestyle="--", marker=".", color="blue", label="Avec Path-Relinking & max_elite = $max_elite")

    legend()
end


function simulation_test(cost, liaisons_contraintes, liaisons_variables, α, nb_iter, max_elite,N_α,sol_opt)
    target = 0.9 * sol_opt
    time_grasp = Vector{Float64}(undef,0)
    time_pr = Vector{Float64}(undef,0)
    time_reactive= Vector{Float64}(undef,0)
    for i in 1:nb_iter
        t1 = grasp(cost,liaisons_contraintes,liaisons_variables,α,100)[4]
        push!(time_grasp,t1)
        t2= grasppr(cost, liaisons_contraintes, liaisons_variables, α, 50, 3,target)[3]
        push!(time_pr,t2)
        t3= reactive_grasp(cost, liaisons_contraintes, liaisons_variables, [0.3,0.6,0.9], 50, 5,target)[5]
        push!(time_reactive,t3)
    end
    time_grasp = sort(time_grasp)
    time_pr = sort(time_pr)
    time_reactive = sort(time_reactive)
    probabilites_grasp = [(i-0.5)/nb_iter for i in 1:nb_iter]
    figure("comparaison", figsize=(6, 6)) # Create a new figure
    title("Time to get 98% of max value ")
    xlabel("temps")
    ylabel("probabilités cumulées")
    xlim(0,5)
    ylim(0, 1)

    plot(time_grasp,probabilites_grasp,marker=".", color="red",label="grasp")
    plot(time_pr,probabilites_grasp,marker=".", color="green",label="path_relinking")
    plot(time_reactive,probabilites_grasp,marker=".", color="blue",label="reactive")
    legend()
end
