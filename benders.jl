using CSV
using CPLEX
using DataFrames
using FacilityLocationProblems
using Gurobi
using JuMP
using HiGHS

include("utils.jl")

global solver = HiGHS  #Gurobi, CPLEX, HiGHS
const GRB_ENV = Gurobi.Env()

function p_median_relaxed_benders(data::PMedianProblem)
    # - Full model relaxed optimal solution 
    z_opt = full_model_relaxed(data)
    
    # - Master Problem 
    master = build_master_problem(data)
    
    # - Subproblem
    sp = build_subproblem_relaxed(data)
    
    optimal = false 
    i = 1 

    # - Initialize variables 
    LB = 0.0
    UB = 0.0

    # DataFrame to store iteration data

    df = DataFrame(iteration=Int[], max_coeff=Float64[],min_coeff=Float64[], constant=Float64[], LB=Float64[], UB=Float64[])
    n_medians = data.medians
    for k in 1:n_medians
        df[!, Symbol("med$k")] = Int[]
    end

    # Saving solutions in a Dict: Iteration -> (x_sol, y_sol)
    dict_solutions = Dict{Int,Tuple{Vector{Int64},Matrix{Float64}}}()

    # dict_solutions = Dict{Int,Vector{Int}}()

    while !optimal
        println("\n -- Iteration: $i -- ") 
        
        # -- Step 1: Solve the master problem  
        x_fixed = solve_master_problem(master) 
        LB = objective_value(master)
        
        # -- Step 2: Update and Solve the subproblem
        set_solution!(data,sp,x_fixed)
        
        coeffs,constant = solve_subproblem(sp,data)
        y_fixed = value.(sp[:y])

        UB = objective_value(sp)
        
        # -- Step 3: Update the master problem
        add_optimality_cut!(master,data,coeffs,constant)
        
        # -- Step 4: Check convergence
        if abs(LB-UB) < 1e-6
            optimal = true
            print("Convergence reached \n")
        end
        
        println("UB: ",UB)
        println("LB: ",LB)
        println("opt:",z_opt)

        # Store solutions 
        dict_solutions[i] = (x_fixed, y_fixed)

        # Store iteration data
        coeffs   = round.(coeffs, digits=2)
        constant = round(constant, digits=2)
        LB       = round(LB, digits=2)
        UB       = round(UB, digits=2)

        
        medians = [ j for j in 1:nc(data) if x_fixed[j] > 0 ]
        medians = sort(medians)
        # print(medians)

        # saving solutions 
        push!(df, vcat(i, maximum(coeffs),minimum(coeffs), constant, LB, UB, medians))

        i += 1
    end
    # Save DataFrame to CSV
    CSV.write("benders_iterations.csv", df, delim = ';')
    return dict_solutions
end

function p_median_relaxed_benders_fixed(data::PMedianProblem, dict_solutions::Dict{Int,Tuple{Vector{Int64},Matrix{Float64}}})
    # - Full model relaxed optimal solution 
    z_opt = full_model_relaxed(data)
    
    # - Master Problem 
    master = build_master_problem_fixed(data)
    
    # - Subproblem
    sp = build_subproblem_relaxed(data)
    
    # - Initialize variables 
    LB = 0.0
    UB = 0.0

    df = DataFrame(iteration=Int[], max_coeff=Float64[],min_coeff=Float64[], constant=Float64[], LB=Float64[], UB=Float64[])
    n_medians = data.medians
    for k in 1:n_medians
        df[!, Symbol("med$k")] = Int[]
    end

    n_iterations = length(dict_solutions)

    for i in 1:n_iterations
        println("\n -- Iteration: $i -- ") 
        
        # -- Get fixed solutions 
        x_fixed = dict_solutions[i][1]
        # y_fixed = dict_solutions[i][2]

        # -- Step 1: Solve the master problem
        set_fixed_master_solution!(data,master,x_fixed)
        
        solve_master_problem(master) 
        LB = objective_value(master)
        
        # -- Step 2: Update and Solve the subproblem
        set_solution!(data,sp,x_fixed)
        
        coeffs,constant = solve_subproblem(sp,data)

        UB = objective_value(sp)
        
        # -- Step 3: Update the master problem
        add_optimality_cut!(master,data,coeffs,constant)
        
        # -- Step 4: Check convergence
        if abs(LB-UB) < 1e-6
            print("Convergence reached \n")
        end
        
        println("UB: ",UB)
        println("LB: ",LB)
        println("opt:",z_opt)

        # Store iteration data
        coeffs   = round.(coeffs, digits=2)
        constant = round(constant, digits=2)
        LB       = round(LB, digits=2)
        UB       = round(UB, digits=2)
        
        medians = [ j for j in 1:nc(data) if x_fixed[j] > 0 ]
        medians = sort(medians)

        # saving solutions 
        push!(df, vcat(i, maximum(coeffs),minimum(coeffs), constant, LB, UB, medians))

    end

    # Save DataFrame to CSV
    CSV.write("benders_iterations_fixed.csv", df, delim = ';')
end

function p_median_mixed_combinatorial_benders(data::PMedianProblem)
    # - Full model relaxed optimal solution 
    z_opt = full_model(data)
    
    # - Master Problem 
    master = build_master_problem(data)
    
    # - Subproblem
    sp = build_subproblem_relaxed(data)
    
    optimal = false 
    i = 1 
    I = 1:nc(data)

    # - Initialize variables 
    LB = 0.0
    UB = 0.0

    π = zeros(nc(data))
    β = zeros(nc(data),nc(data))
    λ = zeros(nc(data))

    integral = false 
    L = 0.0

    while !optimal
        println("\n -- Iteration: $i -- ") 
        
        # --  Solve the master problem  
        x_fixed = solve_master_problem(master) 
        LB = objective_value(master)
        
        # --  Update and Solve the subproblem
        set_solution!(data,sp,x_fixed)
        
        if !integral
            π , β ,λ = solve_subproblem(sp)
            UB = objective_value(sp)
            
            # -- Update the master problem
            add_optimality_cut!(master,data,π,β,λ)
        else 
            optimize!(sp)
            UB = objective_value(sp)
            add_combinatorial_cut!(master,sp,data,x_fixed,L)
        end

        # -- Step 4: Check convergence
        if abs(LB-UB) < 1e-6
            if !integral

                println(" - #################### - ")
                println(">>> Going to Integral! <<<")
                println(" - #################### - ")
                
                for i in I, j in I
                    set_binary(sp[:y][i,j])
                end

                integral = true
            else 
                optimal = true
                print("Convergence reached \n")
            end
        end
        
        println("UB: ",UB)
        println("LB: ",LB)
        println("opt:",z_opt)
        i += 1
    end

end

data = loadPMedianProblem(:pmedcap03)

# dict_solutions = p_median_relaxed_benders(data)

# global solver = Gurobi  #Gurobi, CPLEX, HiGHS

# p_median_relaxed_benders_fixed(data, dict_solutions)



