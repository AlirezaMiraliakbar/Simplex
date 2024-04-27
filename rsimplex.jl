#=
Implementation of revised simplex algorithm including Phase I & II

You should insert only the main problem's matrices and vectors and then the algorithm will take take care of the rest by 
first augmenting the auxilary problem parameters and solving it to give the initial basic feasible solution (BFS) for the 
actual LP problem

author: Alireza Miraliakbar
=#
using LinearAlgebra
function revised_simplex(A, b, c, x_init,indx_B)
    #=
    Ax=b
    -------------------------------------------
    A: the coefficient matrix of the LP problem
    b: the RHS vector of the problem
    c: the cost vector 
    =#

    m,n = size(A)
    # building a tableau
    tableau = zeros(Float64, m+1, n+1)    
    # selecting basic variables 
    x_B = x_init[indx_B]
    # println("x_B = $x_B")
    # selecting basic cost coeffs
    c_B = c[indx_B]
    # println("c_B = $c_B")

    # constructing the initial tableau
    tableau[1,1] = - dot(c_B,x_B)
    c_bar = transpose(c) - transpose(c_B) * A
    # println("c_bar = $c_bar")
    tableau[1,2:n+1] = c_bar
    
    tableau[2:m+1, 1] = transpose(x_B)

    tableau[2:m+1, 2:n+1] = A
    x_opt = 1
    println("-------------------- Initial Tableau -------------------------------")

    display(tableau)
    i = 0
    max_iter = 20
    while i < max_iter
        neg_indx = findall(x -> x < 0, tableau[1,2:n+1])
        # Dantzig's rule
        pivot_indx = neg_indx[1]
        println("pivot_indx = $pivot_indx")
        u_j = tableau[2:m+1,pivot_indx+1]
        println("u_j = $u_j")
        non_zero_indx = findall(x -> x > 0, u_j)
        println("non_zero_indx = $non_zero_indx")
        thetas = tableau[2:m+1,1][non_zero_indx] ./ u_j[non_zero_indx]
        println("thetas = $thetas")
        min_theta = minimum(thetas)
        println("minimum theta = $min_theta")
        l_prime = findall(x -> x == min_theta, thetas)[1]
        l = non_zero_indx[l_prime]
        println("l = $l")
        pivot_element = tableau[l+1,pivot_indx+1]
        println("pivot_element = $pivot_element")
        # making reduced cost of j = 0 
        c_j = tableau[1,2:n+1][pivot_indx]
        # println("c_j = $c_j")
        
        # updating other reduced costs
        coeff_c_bar_j = c_j / pivot_element
        tableau[1,:] += -tableau[l+1,:] * coeff_c_bar_j 
        # make column A_[pivot_indx] into e_l
        
        # println("pivot element = $pivot_element")
        e_l, coeffs = transform_to_unit_vector(u_j, l)
        # println("e_l = $e_l")
        # println("coeffs = $coeffs")
        tableau_prime = apply_coefficients(tableau[2:m+1,1:n+1], coeffs, l)
        tableau_prime[l,pivot_indx+1] = 1
        tableau[2:m+1,1:n+1] = tableau_prime
        i+= 1
        println("-------------------- Iteration $i -------------------------------")

        display(tableau)

        if all(tableau[1,2:n+1] .>= 0.0)
            break
        end
    end
    x_opt = tableau[2:m+1,1]
    return x_opt
end

function transform_to_unit_vector(v, k)
    n = length(v)
    if v[k] == 0.0
        error("Element at index $k is zero, cannot scale.")
    end

    # Initialize the transformation coefficients array.
    coeffs = zeros(n)  # This will store the coefficients used to zero out other entries.

    # First, scale the vector such that the k-th element becomes 1.
    scale_factor = 1 / v[k]
    v .= v .* scale_factor  # Scale the entire vector.

    # Now, zero out all other elements using the k-th element.
    for i = 1:n
        if i != k
            coeffs[i] = -v[i]  # Store the coefficient used for zeroing this element.
            v[i] = 0  # Explicitly zero out the element.
        end
    end

    v[k] = 1  # Ensure the k-th element is exactly 1.
    return v, coeffs
end

function apply_coefficients(matrix, coeffs, row_index)
    n = size(matrix, 1)
    if length(coeffs) != n
        error("Number of coefficients must match the number of rows in the matrix.")
    end
    if row_index < 1 || row_index > n
        error("Invalid row index.")
    end

    # Create a copy of the input matrix to avoid modifying it in place.
    result_matrix = copy(matrix)

    # Apply coefficients to the corresponding row of the matrix.
    for i = 1:n
        if i != row_index
            result_matrix[i, :] .+= coeffs[i] * matrix[row_index, :]
        end
    end
    return result_matrix
end

A = [1.0 2.0 3.0 0.0; -1.0 2.0 6.0 0.0; 0.0 4.0 9.0 0.0; 0.0 0.0 3.0 1.0]
b = [3.0 ;2.0; 5.0; 1.0]
c =[1.0; 1.0; 1.0; 1.0]
m, n = size(A)
# Phase I
A_p1 = hcat(A,I)

# we will provide the initial solution to the simplex 

basic_indx = [5, 6, 7, 8]
x_p1 = [0; 0; 0; 0; 3; 2; 5; 1]
cc = fill(0.0,(m,1))
c_p1 = vcat(cc, c)
x_opt_p1 = revised_simplex(A_p1, b, c_p1, x_p1, basic_indx)