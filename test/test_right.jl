@testset "solver_right" begin

    # #####################################
    # First test (5 possible solutions)
    # ####################################
    A = [1.0 0.0 ; 0.0 1.0]
    B = [-1.0 -6.0 ; 2.0 -9.0]
    C = [0.0 12.0 ; -2.0 14.0]

    # Solve the problem with Tolerance
    tol = 1E-6

    # Solve the problem without Line Search
    X = quadraticsolve_right(A,B,C,tol=tol,niter=100,verbose=true) 

    # Compute Q
    Q = (A*X .+ B)*X .+ C
    
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
    X = quadraticsolve_right(A,B,C,tol=tol,niter=100,verbose=true) 

    # Compute Q
    Q = (A*X .+ B)*X .+ C
    
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
    X = quadraticsolve_right(A,B,C,tol=tol,niter=100,verbose=true) 

    # Compute Q
    Q = (A*X .+ B)*X .+ C
    
    @test norm(Q)<=tol

    
    # ###############################################
    # Fourth test 
    # ###############################################
    A = [17.6  1.28  2.89; 
         1.28 0.824  0.413;
         2.89 0.413  0.725]

    B = [7.66 2.45 2.1;
         0.23 1.04 0.223;
         0.6  0.756 0.658]

    C = [121.0 18.9 15.9 ;
         0.0   2.7  0.145;
         11.9 3.64  15.5]

    # Solve the problem with Tolerance
    tol = 1E-6

    # Solve the problem without Line Search
    X = quadraticsolve_right(A,B,C,tol=tol,niter=100,verbose=true) 

    # Compute Q
    Q = (A*X .+ B)*X .+ C
    
    @test norm(Q)<=tol
  
end
