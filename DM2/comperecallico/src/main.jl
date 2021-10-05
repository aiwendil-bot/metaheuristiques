include("../../../libSPP/librarySPP.jl")
include("pretraitement.jl")
include("codeDM2.jl")
# include("/Users/nicolascompere/GitHub/metaheuristiques/DM1/comperecallico/src/codeDM1.jl")

target = "../../../Data"
C, A = loadSPP(string(target, "/", "pb_100rnd0100.dat"))

liaisons_contraintes = vect_contraintes(A)
liaisons_variables = vect_variables(A)

@elapsed grasp(C, liaisons_contraintes, liaisons_variables, 0.6, 200)
