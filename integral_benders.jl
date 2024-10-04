include("benders.jl")

global solver = Gurobi  #Gurobi, CPLEX, HiGHS

function p_median_new_combinatorial_benders(data::PMedianProblem)
    I = 1:nc(data)

    # - Full model relaxed optimal solution 
    z_opt = full_model(data)
    z_opt_relaxed = full_model_relaxed(data)

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

    df = DataFrame(iteration=Int[], max_coeff=[],min_coeff=[], constant=[], LB=Float64[], UB=Float64[])
    n_medians = data.medians
    for k in 1:n_medians
        df[!, Symbol("med$k")] = Int[]
    end

    # Saving solutions in a Dict: Iteration -> (x_sol, y_sol)
    dict_solutions = Dict{Int,Tuple{Vector{Int64},Matrix{Float64}}}()

    integral = false 
    L = 0.0
    
    while !optimal
        println("\n -- Iteration: $i -- ") 
        
        # -- Step 1: Solve the master problem  
        x_fixed = solve_master_problem(master) 
        LB = objective_value(master)
        
        # -- Step 2: Update and Solve the subproblem
        set_solution!(data,sp,x_fixed)
        
        if !integral
            coeffs,constant = solve_subproblem(sp,data)
            UB = objective_value(sp)
            
            # -- Update the master problem
            add_optimality_cut!(master,data,coeffs,constant)
        else 
            optimize!(sp)
            UB = objective_value(sp)
            add_combinatorial_cut!(master,sp,data,x_fixed,L)
        end

        y_fixed = value.(sp[:y])

        UB = objective_value(sp)
        
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
        println("opt_relaxed:",z_opt_relaxed)

        # Store solutions 
        dict_solutions[i] = (x_fixed, y_fixed)

        # Store iteration data
        LB       = round(LB, digits=2)
        UB       = round(UB, digits=2)

        medians = [ j for j in 1:nc(data) if x_fixed[j] > 0 ]
        medians = sort(medians)

        if !integral
            coeffs   = round.(coeffs, digits=2)
            constant = round.(constant, digits=2)
            push!(df, vcat(i, maximum(coeffs),minimum(coeffs), constant, LB, UB, medians))
        else 
            push!(df, vcat(i, "-", "-", "-", LB, UB, medians))
        end
        i += 1
    end
    # Save DataFrame to CSV
    CSV.write("benders_iterations.csv", df, delim = ';')
    return dict_solutions
end

data = loadPMedianProblem(:pmedcap03)

dict_solutions = p_median_new_combinatorial_benders(data);

