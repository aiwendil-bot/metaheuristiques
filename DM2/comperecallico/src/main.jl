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
show(simulation_test(C, liaisons_contraintes, liaisons_variables, 0.7, 20, 5,20,731))
