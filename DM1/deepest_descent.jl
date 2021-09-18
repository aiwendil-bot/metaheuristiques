include("loadSPP.jl")
include("getfname.jl")
include("greedy.jl")

using LinearAlgebra

cd("$(homedir())/github/metaheuristiques")

fname = "DM1/Data/pb_100rnd0100.dat"
cost, matrix = loadSPP(fname)

#=
function deepest_descent(initial_solution)
    pour au moins 2 voisinages :
    repeat
        choisir x' dans le voisinage de x tq z(x") >= z(x') pour tout x" du voisinage
        if z (x ′) > z (x ) then
        x ←x ′;
        end if
until z (x') <= z (x ), ∀x' du voisinage
Output x ;
end
=#
