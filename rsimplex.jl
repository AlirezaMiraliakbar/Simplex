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
    tableau = fill(0,(m+1,n+1))
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

    while all(tableau[1,2:n+1] .<= 0)
        neg_indx = findall(x -> x < 0, tableau[1,2:n+1])
        # Dantzig's rule
        pivot_indx = neg_indx[1]
        println("pivot_indx = $pivot_indx")
        u_j = A[:,pivot_indx]
        println("u_j = $u_j")
        non_zero_indx = findall(x -> x > 0, u_j)
        println("non_zero_indx = $non_zero_indx")
        thetas = x_B[non_zero_indx] ./ u_j[non_zero_indx]
        println("Theta = $thetas")
        min_theta = minimum(thetas)
        l = findall(x -> x == min_theta, thetas)[1]
        println("l = $l")
        # making reduced cost of j = 0 
        c_j = tableau[1,2:n+1][pivot_indx]
        println("c_j = $c_j")
        tableau[1,1] += min_theta * c_j
        tableau[1,pivot_indx+1] = 0

        # make column A_[pivot_indx] into e_l
        pivot_element = tableau[l+1,pivot_indx+1]
        println("pivot element = $pivot_element")
        
        break
    end

    return x_opt
end

A = [1 2 3 0; -1 2 6 0; 0 4 9 0; 0 0 3 1]
b = [3 ;2; 5; 1]
c =[1; 1; 1; 1]
m, n = size(A)
# Phase I
A_p1 = hcat(A,I)

# we will provide the initial solution to the simplex 

basic_indx = [5, 6, 7, 8]
x_p1 = [0; 0; 0; 0; 3; 2; 5; 1]
cc = fill(0,(m,1))
c_p1 = vcat(cc, c)
x_opt_p1 = revised_simplex(A_p1, b, c_p1, x_p1, basic_indx)