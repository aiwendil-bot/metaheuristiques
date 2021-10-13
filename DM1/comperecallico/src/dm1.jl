#= ATTENTION:
   your own folder is considered as the current working directory
   for running your solver
=#

include("../../../libSPP/librarySPP.jl")

include("codeDM1.jl")
include("solver.jl")
include("pretraitement.jl")

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

        t1 = @elapsed x, zInit = greedy_construction(C, liaisons_contraintes,liaisons_variables)
        t2 = @elapsed xbest, zBest = deepest_descent(C, liaisons_contraintes,liaisons_variables, x, zInit)

        #t1 = @elapsed z = setSPP(C, A)
        #println(fnames[instance]," : ", t1," ", z)

        # Saving results -------------------------------------------------------
        println(io, fnames[instance], " ", zInit, " ", zBest, " ", t1, " ", t2, " ", t1+t2)
    end
    close(io)



end

main()
