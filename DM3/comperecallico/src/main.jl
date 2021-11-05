include("tabou.jl")
include("tabouv2.jl")
include("tabouv3.jl")
include("pretraitement.jl")
include("../../../libSPP/librarySPP.jl")

function main()
    println("Etudiants : COMPERE Nicolas et CALLICO Adrien")

    # Collecting the names of instances to solve located in the folder Data ----
    target = "../../../datatest"
    fnames = getfname(target)

    fres = splitdir(splitdir(pwd())[end-1])[end]
    io = open("../res/"*fres*".res", "w")
    for instance = 1:length(fnames)

        # Load one numerical instance ------------------------------------------
        C, A = loadSPP(string(target,"/",fnames[instance]))
        liaisons_contraintes = vect_contraintes(A)
        liaisons_variables = vect_variables(A)

        zInit = 0 ; zBest = 0 ; t1 =0.0 ; t2 = 0.0

        #un, zInit, probamax, quatre = reactive_grasp(C, liaisons_contraintes, liaisons_variables, vector_α, nb_iter, 20)

        x, zInit = tabouv3(C,liaisons_contraintes,liaisons_variables,10,7)
        meilleurs = [4,372,203,40,184,1004,571,926,122,1141,2236,424]
        #t1 = @elapsed z = setSPP(C, A)
        #println(fnames[instance]," : ", t1," ", z)

        # Saving results -------------------------------------------------------
        println(io, fnames[instance], " ", zInit/meilleurs[instance])
        #show(reactive_experiment(C, liaisons_contraintes, liaisons_variables, vector_α, nb_iter, N_α))
    end
    close(io)



end
