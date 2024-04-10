import numpy as np

def rsimplex(c, A, b):
    # First step is to have the initial basic feasible solution
    m = A.shape[0]
    n = A.shape[1]
    reduced_costs = np.zeros(n)
    # (1) initial BFS by slack variable
    zero_vector = np.zeros([n,1])
    x_bfs_0 = np.vstack([zero_vector, b])
    B = np.identity(m)
    B_inv = np.linalg.inv(B)

    # (2) computing row vector p'
    p_transposed = np.matmul(np.transpose(c), B_inv)

    # (3) computing reduced costs 
    reduced_costs = c - np.dot(p_transposed, A)
    print(reduced_costs)
    
    
    return 1

if __name__ == "__main__":
    A = np.array([[2.43, 2.07], [2.57, 2.44], [7.16, 1.26]])
    # Define vector b
    b = np.array([[2.5447883],
                [6.14664116],
                [3.18367944]])
    # Define vector c
    c = np.array([[4.32838813],
                [9.58224293]])
    print("Matrix A:\n", A)
    print("\nVector b:\n", b)
    print("\nVector c:\n", c)

    results = rsimplex(A,b,c)