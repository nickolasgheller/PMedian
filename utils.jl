

function full_model_relaxed(data::PMedianProblem)::Float64
    p = data.medians
    I = 1:nc(data)
    
    w = data.demands
    d = data.costs
    Q = data.capacity

    if solver == Gurobi
        model = Model(() -> Gurobi.Optimizer(GRB_ENV))
    else
        model = Model(CPLEX.Optimizer)
    end
    
    @variable(model, x[j in I], Bin)
    
    @variable(model,  y[i in I, j in I] >= 0)

    @objective(model, Min, sum(sum(d[i, j]y[i, j] for j in I) for i in I))
    
    @constraint(model, [i in I], sum(y[i, j] for j in I) == 1)
    @constraint(model, [i in I, j in I], y[i, j] <= x[j])
    @constraint(model, sum(x[j] for j in I) == p)
    @constraint(model, [j in I], sum(w[i]y[i, j] for i in I) <= Q)
    
    # print(model)
    set_silent(model)
    optimize!(model)
    return objective_value(model)
end

function full_model(data::PMedianProblem)::Float64
    p = data.medians
    I = 1:nc(data)
    
    w = data.demands
    d = data.costs
    Q = data.capacity

    if solver == Gurobi
        model = Model(() -> Gurobi.Optimizer(GRB_ENV))
    else
        model = Model(CPLEX.Optimizer)
    end
    
    @variable(model, x[j in I], Bin)
    @variable(model, y[i in I, j in I], Bin)
    
    @objective(model, Min, sum(sum(d[i, j]y[i, j] for j in I) for i in I))
    
    @constraint(model, [i in I], sum(y[i, j] for j in I) == 1)
    @constraint(model, [i in I, j in I], y[i, j] <= x[j])
    @constraint(model, sum(x[j] for j in I) == p)
    @constraint(model, [j in I], sum(w[i]y[i, j] for i in I) <= Q)
    
    # print(model)
    set_silent(model)
    optimize!(model)
    return objective_value(model)
end

function build_master_problem(data::PMedianProblem)::Model
    p = data.medians
    n = nc(data)
    
    I = 1:n
    
    if solver == Gurobi
        model = Model(() -> Gurobi.Optimizer(GRB_ENV))
    else
        model = Model(CPLEX.Optimizer)
    end
    
    @variable(model, x[j in I], Bin)
    @variable(model, α >= 0)

    @objective(model, Min, α)
    
    @constraint(model, sum(x[j] for j in I) == p)
    
    # print(model)
    set_silent(model)
    return model
end

function build_master_problem_fixed(data::PMedianProblem)::Model
    p = data.medians
    n = nc(data)
    
    I = 1:n
    
    if solver == Gurobi
        model = Model(() -> Gurobi.Optimizer(GRB_ENV))
    else
        model = Model(CPLEX.Optimizer)
    end
    
    @variable(model, x[j in I], Bin)
    @variable(model, α >= 0)

    @objective(model, Min, α)
    
    @constraint(model, sum(x[j] for j in I) == p)
    
    # @constraint(model, x_fix[j in I], x[j] == 0)

    # print(model)
    set_silent(model)
    return model
end



function build_subproblem_relaxed(data::PMedianProblem)::Model
    I = 1:nc(data)
    
    w = data.demands
    d = data.costs
    Q = data.capacity

    if solver == Gurobi
        model = Model(() -> Gurobi.Optimizer(GRB_ENV))
    else
        model = Model(CPLEX.Optimizer)
    end
    
    @variable(model, y[i in I, j in I] >= 0)
    
    # @variable(model, y[i in I, j in I], Bin)

    @objective(model, Min, sum(sum(d[i, j]y[i, j] for j in I) for i in I))
    
    @constraint(model, bind[i in I, j in I], y[i, j] <= 0)         #β[i,j]
    @constraint(model, customers[i in I], sum(y[i, j] for j in I) == 1)     #π[i] 
    @constraint(model, knapsacks[j in I], sum(w[i]y[i, j] for i in I) <= Q) #λ[j]

    # print(model)
    set_silent(model)
    return model
end

function set_solution!(data::PMedianProblem,subproblem::Model,x_fixed::Vector{Int64})
    I = 1:nc(data)

    for j in I, i in I
        set_normalized_rhs(subproblem[:bind][i,j], x_fixed[j])
    end

end

function set_fixed_master_solution!(data::PMedianProblem,master::Model,x_fixed::Vector{Int64})
    I = 1:nc(data)
    x = master[:x]
    for i in I
        JuMP.fix(x[i], x_fixed[i])
    end
end


function add_optimality_cut!(master::Model,data::PMedianProblem,coeffs::Vector{Float64},constant::Float64)
    I = 1:nc(data)
    # Q = data.capacity
    α = master[:α]
    x = master[:x]

    @constraint(master, α ≥ sum(coeffs[j]*x[j] for j in I) + constant)

    # @constraint(master, α >= sum(π[i] for i in I) + sum(β[i,j]*x[j] for i in I, j in I) + sum(λ[j]*Q for j in I))
end

function add_combinatorial_cut!(master::Model,subproblem::Model,data::PMedianProblem,x_fixed::Vector{Int64},L::Float64)
    I = 1:nc(data)
    z = objective_value(subproblem)
    α = master[:α]
    x = master[:x]
    @constraint(master, α ≥ z - (z-L)*(sum(x[i] for i in I if x_fixed[i] == 0) + sum(1 - x[i] for i in I if x_fixed[i] == 1)))
    # @constraint(master, α ≥ z - (z-L)*(sum(x[i] for i in I if x_fixed[i] == 0)))
    # @constraint(master, α ≥ z - (z-L)*(sum(1 - x[i] for i in I if x_fixed[i] == 1)))
end

function solve_master_problem(master::Model)::Vector{Int64}
    optimize!(master)
    x = value.(master[:x])
    return x
end

function solve_subproblem(subproblem::Model,data::PMedianProblem)::Tuple{Vector{Float64},Float64}
    optimize!(subproblem)
    # π = dual.(subproblem[:customers])
    # β = dual.(subproblem[:bind])
    # λ = dual.(subproblem[:knapsacks])
    I = 1:nc(data)

    coeffs = [ sum(dual.(subproblem[:bind][:, j])) for j in I ]
    constant = sum(dual.(subproblem[:customers]))
    constant += data.capacity * sum(dual.(subproblem[:knapsacks]))

    return coeffs,constant
end

# function set_fixed_subproblem_solution!(data::PMedianProblem,subproblem::Model,x_fixed::Vector{Int64},y_fixed::Matrix{Float64})
#     I = 1:nc(data)
#     y = subproblem[:y]
#     for j in I, i in I
#         # JuMP.fix(y[i,j], y_fixed[i,j], force = true)
#         set_normalized_rhs(subproblem[:bind][i,j], x_fixed[j])
#         # set_normalized_rhs(subproblem[:y_fix][i,j], y_fixed[i,j])
#     end
# end


# function build_subproblem_fixed(data::PMedianProblem)::Model
#     I = 1:nc(data)
    
#     w = data.demands
#     d = data.costs
#     Q = data.capacity

#     if solver == Gurobi
#         model = Model(() -> Gurobi.Optimizer(GRB_ENV))
#     else
#         model = Model(CPLEX.Optimizer)
#     end
    
#     # @variable(model, y[i in I, j in I], Bin)
#     @variable(model, y[i in I, j in I] >= 0)
    
#     @objective(model, Min, sum(sum(d[i, j]y[i, j] for j in I) for i in I))
    
#     @constraint(model, bind[i in I, j in I], y[i, j] <= 0)         #β[i,j]
#     @constraint(model, customers[i in I], sum(y[i, j] for j in I) == 1)     #π[i] 
#     @constraint(model, knapsacks[j in I], sum(w[i]y[i, j] for i in I) <= Q) #λ[j]

#     # @constraint(model, y_fix[i in I, j in I], y[i,j] == 0)
#     # print(model)
#     set_silent(model)
#     return model
# end