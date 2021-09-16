include("loadSPP.jl")
include("getfname.jl")

cd("$(homedir())/github/metaheuristiques")
#test ouvrir un fichier d'instance
fname = "DM1/Data/didactic.dat"
cost, matrix = loadSPP(fname)
println(cost)
println(matrix)
