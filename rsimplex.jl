#=
Implementation of revised simplex algorithm including Phase I & II

You should insert only the main problem's matrices and vectors and then the algorithm will take take care of the rest by 
first augmenting the auxilary problem parameters and solving it to give the initial basic feasible solution (BFS) for the 
actual LP problem

author: Alireza Miraliakbar
=#
using LinearAlgebra
function revised_simplex(A::Matrix, b::Vector, c::Vector, phase::String)
    #=
    Ax=b
    -------------------------------------------
    A: the coefficient matrix of the LP problem
    b: the RHS vector of the problem
    c: the cost vector 
    =#

    m,n = size(A)
    if phase == "phase I"
        B_init = Matrix{Float64}(I, m,m)
        x_B = b
        c_bar = fill(0,(n,1))
        c_B = fill(1,(n,1))
        c_bar_nonbasic = - transpose(c_B) * A
        c_bar_basic = fill(1,(n,1))
        tableau = fill(0,(m,n+1))
        tableau[0,1:]
        while all(x -> x >= 0, c_bar)
            
        end
    end
    return c_bar
end
# test

# b = []
# # Phase I
A= [1 2 3 0; -1 2 6 0; 0 4 9 0; 0 0 3 1]
b = [2 ;3]
c =[1; 1; 1]
x_opt_p1 = revised_simplex(A, b, c, "phase I")
# # Phase II
# A= [1 2 3 0; -1 2 6 0; 0 4 9 0; 0 0 3 1]
# m, n = size(A)

# c_p2 = []
# x_opt_p2 = revised_simplex(A_p2, b, c_p2)