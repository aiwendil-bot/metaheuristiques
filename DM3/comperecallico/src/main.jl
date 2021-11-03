include("tabou.jl")
include("tabouv2.jl")
include("tabouv3.jl")
include("pretraitement.jl")
include("../../../libSPP/librarySPP.jl")

target = "../../../Data"
C, A = loadSPP(string(target, "/", "pb_100rnd0100.dat"))

liaisons_contraintes = vect_contraintes(A)
liaisons_variables = vect_variables(A)

tabouv3(C,liaisons_contraintes,liaisons_variables,10,50)
