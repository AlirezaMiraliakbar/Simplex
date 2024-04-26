#=
Implementation of revised simplex algorithm including Phase I & II

You should insert only the main problem's matrices and vectors and then the algorithm will take take care of the rest by 
first augmenting the auxilary problem parameters and solving it to give the initial basic feasible solution (BFS) for the 
actual LP problem

author: Alireza Miraliakbar
=#
using LinearAlgebra
function revised_simplex(A, b, c)
    #=
    Ax=b
    -------------------------------------------
    A: the coefficient matrix of the LP problem
    b: the RHS vector of the problem
    c: the cost vector 
    =#
    m,n = size(A)
    B_init = Matrix{Float64}(I, m,m)
    x_opt = 1
    return x_opt
end
# test
x_test = revised_simplex([1 2 3;4 5 5], [2 ;3], [1 1 1])

# b = []
# # Phase I
# A_p1 = []
# c_p1 = []
# x_opt_p1 = revised_simplex(A_p1, b, c_p1)
# # Phase II
# A_p2 = []
# c_p2 = []
# x_opt_p2 = revised_simplex(A_p2, b, c_p2)