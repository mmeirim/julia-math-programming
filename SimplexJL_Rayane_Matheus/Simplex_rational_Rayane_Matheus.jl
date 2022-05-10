using LinearAlgebra, Combinatorics, Printf

mutable struct SimplexTableau
    z_c::Array{Rational{Int64}} # z_j - c_j
    Y::Array{Rational{Int64}} # inv(B) * A
    x_B::Array{Rational{Int64}} # inv(B) * b
    obj::Rational{Int64}        # c_B * x_B
    b_idx::Array{Int64}   # indices for basic variables x_B
end

function is_nonnegative(x::Vector)
    return length(x[x.<0]) == 0
end

function is_nonnegativearray(z::Array)
    return length(z[z.<0]) == 0
end

function print_tableau(t::SimplexTableau)
    m, n = size(t.Y)

    hline0 = repeat("-", 6)
    hline1 = repeat("-", 7*n)
    hline2 = repeat("-", 7)
    hline = join([hline0, "+", hline1, "+", hline2])

    println(hline)

    @printf("%6s|", "")
    for j in 1:length(t.z_c)
      @printf("%6s ", string(t.z_c[j]))
    end
    @printf("| %6s\n", string(t.obj))

    println(hline)

    for i in 1:m
      @printf("x[%2d] |", t.b_idx[i])
      for j in 1:n
        @printf("%6s ", string(t.Y[i,j]))
      end
      @printf("| %6s\n", string(t.x_B[i]))
    end

    println(hline)
  end

function pivoting!(t::SimplexTableau)
    m, n = size(t.Y)

    entering, exiting = pivot_point(t)
    println("Pivoting: entering = x_$entering, exiting = x_$(t.b_idx[exiting])")

    # Pivoting: exiting-row, entering-column
    # updating exiting-row
    coef = t.Y[exiting, entering]
    t.Y[exiting, :] /= coef
    t.x_B[exiting] /= coef

    # updating other rows of Y
    for i in setdiff(1:m, exiting)
        coef = t.Y[i, entering]
        t.Y[i, :] -= coef * t.Y[exiting, :]
        t.x_B[i] -= coef * t.x_B[exiting]
    end

    # updating the row for the reduced costs
    coef = t.z_c[entering]
    t.z_c -= coef * t.Y[exiting, :]'
    t.obj -= coef * t.x_B[exiting]

    # Updating b_idx
    t.b_idx[findfirst(t.b_idx .== t.b_idx[exiting])] = entering
end

function pivot_point(t::SimplexTableau)
    # Finding the entering variable index
    entering = findfirst(t.z_c .< 0)[2]
    if entering == 0
        error("Optimal")
    end

    # min ratio test / finding the exiting variable index
    pos_idx = findall(t.Y[:, entering] .> 0)
    if length(pos_idx) == 0
        error("Unbounded")
    end
    exiting = pos_idx[argmin(t.x_B[pos_idx] ./ t.Y[pos_idx, entering])]

    return entering, exiting
end

function initialize(c, A, b, b_idx)
    c = Array{Rational{Int64}}(c)
    A = Array{Rational{Int64}}(A)
    b = Array{Rational{Int64}}(b)

    m, n = size(A)

    # Finding an initial BFS
    # b_idx, x_B, B = initial_BFS(A, b)
    B = A[:, b_idx]
    x_B = inv(B) * b
    Y = inv(B) * A
    c_B = c[b_idx]
    obj = dot(c_B, x_B)

    # z_c is a row vector
    z_c = zeros(1, n)
    n_idx = setdiff(1:n, b_idx)
    z_c[n_idx] = c_B' * inv(B) * A[:, n_idx] - c[n_idx]'
    return SimplexTableau(z_c, Y, x_B, obj, b_idx)
end

function is_optimal(t::SimplexTableau)
    return is_nonnegativearray(t.z_c)
end

function is_unbounded(t::SimplexTableau)
    pos_idx = findall(t.Y[:, entering] .> 0)
    return length(pos_idx) == 0
end

function simplex(A::Matrix{Rational{Int64}}, b::Vector{Rational{Int64}}, c::Vector{Rational{Int64}})
    x = zeros(Rational{Int64}, length(c))

    # Adding slack variables
    A = [A I]
    b_idx = collect(length(c)+1 :size(A,2)) # calculating number of slack variables add and adding their idx to the base
    c = [c;zeros(size(A,1))]
    # SIMPLEX!

    tableau = initialize(c, A, b, b_idx)
    print_tableau(tableau)
    while !is_optimal(tableau)
        pivoting!(tableau)
        print_tableau(tableau)
    end

    if is_optimal(tableau)  # Solução ótima
        opt_x = zeros(Rational{Int64}, length(c))
        opt_x[tableau.b_idx] = tableau.x_B
        println("OPTIMAL SOLUTION: Obj.Value = ",tableau.obj, ", and x = ", opt_x )
        return tableau.obj, x
    elseif is_unbounded(tableau) # Solução ilimitada
        println("UNBOUNDED SOLUTION")
        return Inf64, nothing
    else # Solução inviável
        println("INFEASIBLE SOLUTION")
        return nothing, nothing
    end

end

# Problema do slide do Simplex
# A = [
#     2//1 1//1
#     1//1 2//1
# ]
# b = [ 4//1, 4//1 ]
# c = [ 4//1, 3//1 ]

A = [
    2//1 1//1 1//1 0//1;
    1//1 3//1 0//1 1//1;
]
b = [ 6//1, 9//1]
c = [ 3//1, 1//1, 0//1, 0//1]

A = [
    0//1 1//1
    1//1 1//1
    1//1 0//1
    5//1 1//1
]
b = [ 4//1, 6//1, 3//1, 18//1 ]
c = [ 1//1, 3//1]

simplex(A, b, c)