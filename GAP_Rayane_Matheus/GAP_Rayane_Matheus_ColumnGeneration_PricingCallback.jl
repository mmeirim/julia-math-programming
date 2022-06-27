using AssignmentProblems, JuMP, Gurobi, DataFrames, CSV, BlockDecomposition, Coluna, MathOptInterface, Knapsacks
const BD = BlockDecomposition;
const MOI = MathOptInterface;

global data = AssignmentProblem;
global model = BlockModel;

function criaGAP(data::AssignmentProblem)
    coluna = optimizer_with_attributes(
        Coluna.Optimizer,
        "params" => Coluna.Params(
            solver = Coluna.Algorithm.TreeSearchAlgorithm(maxnumnodes=3) # default branch-cut-and-price
        ),
        "default_optimizer" => optimizer_with_attributes(Gurobi.Optimizer) # Gurobi for the master & the subproblems
    );
    
    model = BlockModel(coluna)

    @axis(M_axis, 1:na(data));
    #resolve o modelo
    @variable(model, x[m in M_axis, j in 1:nj(data)], Bin);
    @constraint(model, cov[j in 1:nj(data)], sum(x[m, j] for m in M_axis) >= 1);
    # @constraint(model, knp[m in M_axis], sum(data.consumptions[m, j] * x[m, j] for j in 1:nj(data)) <= data.capacities[m]);
    @objective(model, Min, sum(data.costs[m, j] * x[m, j] for m in M_axis, j in 1:nj(data)));
    
    @dantzig_wolfe_decomposition(model, decomposition, M_axis)
    subproblems = BD.getsubproblems(decomposition);
    BD.specify!.(subproblems, lower_multiplicity = 0, solver = my_pricing_callback);
        
    return model
end


function resolveGAP(model::Model)
    optimize!(model)
    if has_values(model)
        lb = objective_bound(model)
        ub = objective_value(model)
        tempo = solve_time(model)
    else 
        lb = -Inf
        ub = Inf
        tempo = solve_time(model)
    end
    return lb, ub, tempo, ""
end

# function solve_knapsack(cost, weight, capacity)
#     sp_model = Model(Gurobi.Optimizer)
#     # set_silent(sp_model)
#     items = 1:length(weight)
#     @variable(sp_model, x[i in items], Bin)
#     @constraint(sp_model, weight' * x <= capacity)
#     @objective(sp_model, Min, cost' * x)
#     optimize!(sp_model)
#     x_val = value.(x)
#     return objective_value(sp_model), filter(i -> x_val[i] ≈ 1, collect(items))
# end

function my_pricing_callback(cbdata)
    # Retrieve the index of the subproblem (it will be one of the values in M_axis)
    cur_machine = BD.callback_spid(cbdata, model)

    # Uncomment to see that the pricing callback is called.
    # println("Pricing callback for machine $(cur_machine).")

    # Retrieve reduced costs of subproblem variables
    red_costs = [BD.callback_reduced_cost(cbdata, model[:x][cur_machine, j]) for j in 1:nj(data)]

    profits = deepcopy(red_costs)
    for j in 1:nj(data)
        profits[j] = trunc(-profits[j] * 10^5)
        if profits[j] < 0
            profits[j] = 0
        end
    end

    profits = convert.(Int64, profits)

    # Run the knapsack algorithm
    optimal, jobs_assigned_to_cur_machine = solveKnapsack(Knapsack(data.capacities[cur_machine], data.consumptions[cur_machine, :],profits), :BinaryModel;optimizer = optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag" =>0))

    # Create the solution (send only variables with non-zero values)
    sol_vars = [model[:x][cur_machine, j] for j in jobs_assigned_to_cur_machine]
    sol_vals = [1.0 for _ in jobs_assigned_to_cur_machine]
    sol_cost = sum(red_costs[j] for j in jobs_assigned_to_cur_machine)
    
    # Submit the solution to the subproblem to Coluna
    MOI.submit(model, BD.PricingSolution(cbdata), sol_cost, sol_vars, sol_vals)

    # Submit the dual bound to the solution of the subproblem
    # This bound is used to compute the contribution of the subproblem to the lagrangian
    # bound in column generation.
    MOI.submit(model, BD.PricingDualBound(cbdata), sol_cost) # optimal solution
    return
end

instancias = [
    "a05100",
    "a05200",
    "a10100",
    "a10200",
    "a20100",
    "a20200",
    "b05100",
    "b05200",
    "b10100",
    "b10200",
    "b20100",
    "b20200",
    "c05100",
    "c0515_1",
    "c0515_2",
    "c0515_3",
    "c0515_4",
    "c0515_5",
    "c05200",
    "c0520_1",
    "c0520_2",
    "c0520_3",
    "c0520_4",
    "c0520_5",
    "c0525_1",
    "c0525_2",
    "c0525_3",
    "c0525_4",
    "c0525_5",
    "c0530_1",
    "c0530_2",
    "c0530_3",
    "c0530_4",
    "c0530_5",
    "c0824_1",
    "c0824_2",
    "c0824_3",
    "c0824_4",
    "c0824_5",
    "c0832_1",
    "c0832_2",
    "c0832_3",
    "c0832_4",
    "c0832_5",
    "c0840_1",
    "c0840_2",
    "c0840_3",
    "c0840_4",
    "c0840_5",
    "c0848_1",
    "c0848_2",
    "c0848_3",
    "c0848_4",
    "c0848_5",
    "c10100",
    "c10200",
    "c1030_1",
    "c1030_2",
    "c1030_3",
    "c1030_4",
    "c1030_5",
    "c10400",
    "c1040_1",
    "c1040_2",
    "c1040_3",
    "c1040_4",
    "c1040_5",
    "c1050_1",
    "c1050_2",
    "c1050_3",
    "c1050_4",
    "c1050_5",
    "c1060_1",
    "c1060_2",
    "c1060_3",
    "c1060_4",
    "c1060_5",
    "c15900",
    "c20100",
    "c201600",
    "c20200",
    "c20400",
    "c30900",
    "c401600",
    "c40400",
    "c60900",
    "c801600",
    "d05100",
    "d05200",
    "d10100",
    "d10200",
    "d10400",
    "d15900",
    "d20100",
    "d201600",
    "d20200",
    "d20400",
    "d30900",
    "d401600",
    "d40400",
    "d60900",
    "d801600",
    "e05100",
    "e05200",
    "e10100",
    "e10200",
    "e10400",
    "e15900",
    "e20100",
    "e201600",
    "e20200",
    "e20400",
    "e30900",
    "e401600",
    "e40400",
    "e60900",
    "e801600"
]

results = DataFrame(
   Instance = String[], 
   Model_LB = Float64[],
   Model_UB = Float64[],
   Literature_LB = Float64[],
   Literature_UB = Float64[],
   Solve_time = Float64[],
   Error_Cause = String[]
)

results_line = copy(results);


for instancia in instancias
    println("Starting to solve $instancia")
    solved = Tuple
    try
        global data = loadAssignmentProblem(Symbol(instancia))
        global model = criaGAP(data)

        solved = resolveGAP(model)
    catch e
        solved = (-1, -1, -1., e.msg)
    end

    println("Solved $instancia")
    
    push!(results_line, (instancia, solved[1], solved[2], data.lb, data.ub, solved[3], solved[4]))
    push!(results, (instancia, solved[1], solved[2], data.lb, data.ub, solved[3], solved[4]))
    
    println("$instancia - Bounds from solution: ($(solved[1]), $(solved[2])) | Bounds from literature: ($(data.lb), $(data.ub)) | Time: $(solved[3]) | Error: $(solved[4])") 
    println("\n",results,"\n")
    
    # CSV.write(string("GAP_Rayane_Matheus/","resultsColumnPricingCallback.csv"), results_line, decimal=',', delim=';', append=true)
    empty!(results_line)
end
println(results)
# CSV.write(string("GAP_Rayane_Matheus/","resultsColumnPricingCallbackEnd.csv"), results, decimal=',', delim=';')