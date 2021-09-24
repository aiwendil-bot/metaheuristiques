include("getfname.jl")
include("loadSPP.jl")

using JuMP, GLPK

fname = "DM1/Data/pb_100rnd0100.dat"
cost, matrix = loadSPP(fname)

function model_solve_set_packing(solverSelected::DataType, c::Vector{Int64}, t::Matrix{Int64})
    m::Model = Model(solverSelected)
    nbVar::Int64 = size(t, 2)
    nbConst::Int64 = size(t, 1)
    
    # Variable, objectif et contrainte
    @variable(m, x[1:nbVar] >= 0, binary = true)
    
    @objective(m, Max, sum(c[i]x[i] for i in 1:nbVar))

    @constraint(m, disjoint[j=1:nbConst], sum(t[j,i]x[i] for i in 1:nbVar) <= 1)

    # Résolution
    optimize!(m)

    status = termination_status(m)

    if status == MOI.OPTIMAL
        println("Résolution optimal")
        println("z = ", objective_value(m))
        println("x = ", value.(x))
    elseif status == MOI.INFEASIBLE
        println("Problème impossible")
    elseif status == MOI.INFEASIBLE_OR_UNBOUNDED
        println("Problème non borné")
    end
end

# Appel de la fonction de résolution

function solve_set_packing()
    model_solve_set_packing(GLPK.Optimizer, cost, matrix)
end
