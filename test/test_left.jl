@testset "solver_left" begin

    # #####################################
    # First test (5 possible solutions)
    # ####################################
    A = [1.0 0.0 ; 0.0 1.0]
    B = [-1.0 -6.0 ; 2.0 -9.0]
    C = [0.0 12.0 ; -2.0 14.0]

    # Solve the problem with Tolerance
    tol = 1E-6

    # Solve the problem without Line Search
    X = quadraticsolve_left(A,B,C,LS=false,tol=tol,niter=100,verbose=true) 

    # Compute Q
    Q = X*(X*A .+ B) .+ C
    
    @test norm(Q)<=tol

    # ###############################################
    # Second test (two real and ∞ complex solutions)
    # ###############################################
    A = [1.0 0.0 ; 0.0 1.0]
    B = [-1.0 -1.0 ; 1.0 -1.0]
    C = [0.0 1.0 ; -1.0 0.0]

    # Solve the problem with Tolerance
    tol = 1E-6

    # Solve the problem without Line Search
    X = quadraticsolve_right(A,B,C,LS=false,tol=tol,niter=100,verbose=true) 

    # Compute Q
    Q = X*(X*A .+ B) .+ C
    
    @test norm(Q)<=tol

    # ###############################################
    # Third test (two real and ∞ complex solutions)
    # ###############################################
    A = [1.0 0.0 ; 0.0 1.0]
    B = [1.0 0.0 ; 0.0 1.0]
    C = [-8.0 -12.0 ; -18.0 -26.0]

    # Solve the problem with Tolerance
    tol = 1E-6

    # Solve the problem without Line Search
    X = quadraticsolve_right(A,B,C,LS=false,tol=tol,niter=100,verbose=true) 

    # Compute Q
    Q = X*(X*A .+ B) .+ C
    
    @test norm(Q)<=tol

    
  
end
