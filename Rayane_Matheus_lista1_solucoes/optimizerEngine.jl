using JuMP, GLPK

function runOptimization(oper, nvar, nrest, indices, sinal)
    model = Model(GLPK.Optimizer)

    @variables(model,
    begin
        x[1:nvar] ≥ 0
    end)
    
    for i in 1:nrest
        if sinal[i] == "≤"
            @constraint(model, sum(indices[i,j]*x[j] for j in 1:nvar) ≤ indices[i,end])
        elseif sinal[i] == "≥"
            @constraint(model, sum(indices[i,j]*x[j] for j in 1:nvar) ≥ indices[i,end])   
        end
    end
    
    if oper == "Max"
        @objective(model, Max, sum(x[i]*(indices[end,i]) for i = 1:nvar))
    elseif oper == "Min"
        @objective(model, Min, sum(x[i]*(indices[end,i]) for i = 1:nvar))
    end
    print(model)
    optimize!(model)
    println()
    println(raw_status(model))
    println("z = ", round(objective_value(model),digits=3))
    for j in 1:nvar
        println("x[$j] = ", round(value(x[j]),digits=3))
    end 
end