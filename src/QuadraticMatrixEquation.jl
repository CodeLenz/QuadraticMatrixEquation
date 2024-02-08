module QuadraticMatrixEquation
   
    using LinearAlgebra
    using MatrixEquations
    
    # Load the subrotines
    include("solver_right.jl")
    include("solver_left.jl")

    # Export the main method
    export quadraticsolve_right, quadraticsolve_left

end
