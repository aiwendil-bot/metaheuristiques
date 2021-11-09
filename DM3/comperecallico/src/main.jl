include("tabou.jl")
include("tabouv2.jl")
include("tabouv3.jl")
include("pretraitement.jl")
include("../../../libSPP/librarySPP.jl")

function main()
    println("Etudiants : COMPERE Nicolas et CALLICO Adrien")

    # Collecting the names of instances to solve located in the folder Data ----
    target = "../../../Data"
    fnames = getfname(target)

    fres = splitdir(splitdir(pwd())[end-1])[end]
    io = open("../res/"*fres*".res", "w")
    for instance = 1:length(fnames)

        # Load one numerical instance ------------------------------------------
        C, A = loadSPP(string(target,"/",fnames[instance]))
        liaisons_contraintes = vect_contraintes(A)
        liaisons_variables = vect_variables(A)
        println(greedy_randomized_construction(C, liaisons_contraintes,liaisons_variables, 0.1))
        zInit = 0 ; zBest = 0 ; t1 =0.0 ; t2 = 0.0

        #un, zInit, probamax, quatre = reactive_grasp(C, liaisons_contraintes, liaisons_variables, vector_α, nb_iter, 20)

        x, zInit = tabou(C,liaisons_contraintes,liaisons_variables,10,7)
        #t1 = @elapsed z = setSPP(C, A)
        #println(fnames[instance]," : ", t1," ", z)

        # Saving results -------------------------------------------------------
        println(io, fnames[instance], " ", zInit)
        #show(reactive_experiment(C, liaisons_contraintes, liaisons_variables, vector_α, nb_iter, N_α))
    end
    close(io)



end
