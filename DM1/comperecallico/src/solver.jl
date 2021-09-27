using JuMP, GLPK

function setSPP(C, A)
  m, n = size(A)
  #ip = Model( with_optimizer(GLPK.Optimizer) )
  model = Model(GLPK.Optimizer)
  set_time_limit_sec(model, 300.0)
  @variable(model, x[1:n], Bin)
  @objective(model, Max, dot(C, x))
  @constraint(model , cte[i=1:m], sum(A[i,j] * x[j] for j=1:n) <= 1)
  optimize!(model)
  return objective_value(model)
end
