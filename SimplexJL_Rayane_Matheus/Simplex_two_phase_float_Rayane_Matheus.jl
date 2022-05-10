using LinearAlgebra, Combinatorics, Printf

mutable struct SimplexTableau
    z_c::Array{Float64} # z_j - c_j
    Y::Array{Float64} # inv(B) * A
    x_B::Array{Float64} # inv(B) * b
    obj::Float64        # c_B * x_B
    b_idx::Array{Int64}   # indices for basic variables x_B
end

function is_nonnegative(x::Vector)
    return length(x[x.<0]) == 0
end

function is_nonnegativearray(z::Array)
    return length(z[z.<0]) == 0
end

function is_nonpositive(z::Array)
    return length( z[ z .> 0] ) == 0
end

function initial_BFS(A, b, b_idx, art_idx)
    m, n = size(A)

    B = A[:, b_idx]
    x_B = inv(B) * b

    Y = inv(B) * A
    obj = 0
  
    # z_c is a row vector
    z_c = ones(1,n)
    n_idx = setdiff(1:n, art_idx)
    z_c[n_idx] .= 0
  
    st = SimplexTableau(z_c, Y, x_B, obj, b_idx)   
    print_tableau(st)
    st = adjust_initial_tableau(st, art_idx)

    while !is_optimal(st)
        pivoting!(st)
        print_tableau(st)
    end
    if st.obj < 0
        error("Infeasible")
    end
    return st 
    
end

function adjust_initial_tableau(st::SimplexTableau, art_idx)
    mm, nn = size(st.Y)
    for (aa,gg) in enumerate(art_idx)
        coef = st.Y[aa,gg]
        st.Y[aa, :] /= coef
        st.x_B[aa] /= coef 
        for i in setdiff(1:mm, aa)
            coef = st.Y[i, gg]
            st.Y[i, :] -= coef * st.Y[aa, :]
            st.x_B[i] -= coef * st.x_B[aa]
        end
        # updating the row for the reduced costs
        coef = st.z_c[gg]
        st.z_c -= coef * st.Y[aa, :]'
        st.obj -= coef * st.x_B[aa]
    end
    print_tableau(st)
    return st
end

function print_tableau(t::SimplexTableau)
    m, n = size(t.Y)

    hline0 = repeat("-", 6)
    hline1 = repeat("-", 7 * n)
    hline2 = repeat("-", 7)
    hline = join([hline0, "+", hline1, "+", hline2])

    println(hline)

    @printf("%6s|", "")
    for j in 1:length(t.z_c)
        @printf("%6.2f ", t.z_c[j])
    end
    @printf("| %6.2f\n", t.obj)

    println(hline)

    for i in 1:m
        @printf("x[%2d] |", t.b_idx[i])
        for j in 1:n
            @printf("%6.2f ", t.Y[i, j])
        end
        @printf("| %6.2f\n", t.x_B[i])
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

function initialize(c, A, b, b_idx, art_idx)
    c = Array{Float64}(c)
    A = Array{Float64}(A)
    b = Array{Float64}(b)
  
    m, n = size(A)
  
    # Finding an initial BFS
    println("====================== FIRST PHASE ======================")
    st = initial_BFS(A,b, b_idx, art_idx)
    
    println("====================== SECOND PHASE ======================")
    rm_idx = findall(st.z_c[:] .== 1.0)
    st.z_c = st.z_c[:,setdiff(1:end, rm_idx)]
    st.Y = st.Y[setdiff(1:end, rm_idx), setdiff(1:end, rm_idx)]

    for i in 1: length(st.z_c)
        st.z_c[i] = -c[i]
    end
    
    print_tableau(st)
    adjust_initial_tableau(st,st.b_idx)
    return st
end

function is_optimal(t::SimplexTableau)
    return is_nonnegativearray(t.z_c)
end

function is_unbounded(t::SimplexTableau)
    pos_idx = findall(t.Y[:, entering] .> 0)
    return length(pos_idx) == 0
end

function simplex(A::Matrix{Float64}, b::Vector{Float64}, c::Vector{Float64}, s::Vector{String})
    x = zeros(Float64, length(c))

    # Adding slack, excess and artificial variables
    art_idx = Vector{Int64}()
    slk_idx = Vector{Int64}()
    for s_idx in 1:length(s)
        slk = zeros(Float64, length(s))
        exc = zeros(Float64, length(s))
        art = zeros(Float64, length(s))
        if s[s_idx] == "="
            art[s_idx] = 1
            A = [A art]
            append!(art_idx,size(A,2))
        elseif s[s_idx] == "≤"
            slk[s_idx] = 1
            A = [A slk]
            append!(slk_idx,size(A,2))
        elseif s[s_idx] == "≥"
            exc[s_idx] = -1
            art[s_idx] = 1
            A = [A exc art]
            append!(art_idx,size(A,2))
        end
    end
    
    if length(art_idx) == length(c)
        b_idx = art_idx
    elseif length(art_idx) < length(c)
        diff = length(c) - length(art_idx)
        b_idx = vcat(art_idx,slk_idx[1:diff])
    end
    
    c = [c;zeros(size(A,1))]
    
    # SIMPLEX!

    tableau = initialize(c, A, b, b_idx, art_idx)
    print_tableau(tableau)
    while !is_optimal(tableau)
        pivoting!(tableau)
        print_tableau(tableau)
    end

    if is_optimal(tableau)  # Solução ótima
        opt_x = zeros(length(c))
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

# Problema do slide do Simplex Duas Fases
A = [
    1.0  1.0 1.0
    1.0 -1.0 0.0
    2.0  3.0 1.0
]
b = [ 10.0, 1.0, 20.0 ]
c = [ 4.0, 5.0, -3.0 ]
s = [ "=", "≥", "≤" ]

# A = [
#     0.0 1.0
#     2.0 3.0
# ]
# b = [1.0, 3.0]
# c = [4.0, 2.0]
# s = ["=", "≥"]

simplex(A, b, c, s)