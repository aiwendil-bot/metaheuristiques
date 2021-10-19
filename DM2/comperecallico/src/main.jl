include("../../../libSPP/librarySPP.jl")
include("pretraitement.jl")
include("codeDM2.jl")
include("reactive_grasp.jl")
include("path_relinking.jl")
include("experiment.jl")
# include("/Users/nicolascompere/GitHub/metaheuristiques/DM1/comperecallico/src/codeDM1.jl")

#target = "Data"
target = "../../../Data"
C, A = loadSPP(string(target, "/", "pb_200rnd0300.dat"))

liaisons_contraintes = vect_contraintes(A)
liaisons_variables = vect_variables(A)

# @elapsed println(grasp(C, liaisons_contraintes, liaisons_variables, 0.6, 200))
# @elapsed println(reactive_grasp(C, liaisons_contraintes, liaisons_variables, [0.2, 0.4, 0.6, 0.8], 1000, 20))
# @elapsed println(grasppr(C, liaisons_contraintes, liaisons_variables, 0.6, 100, 4))
#simulation_Reactive()
# reactive_experiment(C, liaisons_contraintes, liaisons_variables, [0.2,0.4,0.6,0.8], 500, 20)
#comparaisons(C, liaisons_contraintes, liaisons_variables, 0.6, 300, 5)


function main()
    println("Etudiants : COMPERE Nicolas et CALLICO Adrien")

    # Collecting the names of instances to solve located in the folder Data ----
    target = "../../../datatest"
    fnames = getfname(target)#[8,11]

    fres = splitdir(splitdir(pwd())[end-1])[end]
    io = open("../res/"*fres*".res", "w")
    for instance = 1:length(fnames)

        # Load one numerical instance ------------------------------------------
        C, A = loadSPP(string(target,"/",fnames[instance]))
        liaisons_contraintes = vect_contraintes(A)
        liaisons_variables = vect_variables(A)

        zInit = 0 ; zBest = 0 ; t1 =0.0 ; t2 = 0.0 ; α = 0.7 ; nb_iter = 200 ; max_elite = 10
        vector_α = [0.2,0.35,0.5,0.65,0.8,0.95]

        #un, zInit, probamax, quatre = reactive_grasp(C, liaisons_contraintes, liaisons_variables, vector_α, nb_iter, 20)
        zInit= grasppr(C, liaisons_contraintes, liaisons_variables, α, nb_iter, max_elite)[2]
        meilleurs = [4,372,203,40,184,1004,571,926,122,1141,2236,424]
        meilleurs_nouveaux = [926,2236]
        #t1 = @elapsed z = setSPP(C, A)
        #println(fnames[instance]," : ", t1," ", z)

        # Saving results -------------------------------------------------------
        println(io, fnames[instance], " ", zInit/meilleurs[instance])
    end
    close(io)



end

main()
