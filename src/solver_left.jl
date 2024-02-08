
# Main function 
#
# Solve the quadratic matrix equation
#
#    X²A + XB + C = 0
#
# where A, B and C are (possible complex) n × n matrices
# and 0 is a null n × n matrix. X is the left solvent.
#
# The Newton-Raphson approach discussed in 
#
# "Solving a quadratic matrix equation by Newton's method 
#  with exact line searchs" by Nicholas J. Higham and
#  Hyun-Min Kim, 2001, SIAM J. Matrix. Anal. Appl. Vol. 23
#  No. 2, pp. 303-316
#
# is used. Contrary to the reference paper, we are using the 
# NonlinearSolve Julia library to solve the line search 
# and MatrixEquation to solve the generalized Sylvester equation
# in each iteration of the NR method.
#
#
# It is important to stress that this problem can be solved by
# calling quadraticsolve_right with A', B' and C' also 
#
#
function quadraticsolve_left(A::AbstractMatrix{T},B::AbstractMatrix{T},C::AbstractMatrix{T};
                             tol=1E-6,niter=100,verbose=false) where T

   #  Problem size
   n = size(A,1)

   # Compute the norms of the input matrices
   nA = norm(A)
   nB = norm(B)
   nC = norm(C)

   # We start by choosing a proper initial 
   # value for solution X
   X =  Matrix((nB + sqrt(nB^2 + 4*nA*nC)/(2*nA))*I(n)*one(Complex{T}))

   # Pre-allocate some matrices used inside the main loop
   XAB = similar(X)
   Q   = similar(X)
   E   = similar(X)

   # Flag for convergence
   cflag = false

   # Intial value for the line-search step
   t = 1.0

   # Main loop
   for iter=1:niter

     # Compute XA - B just once per iteration
     XAB .= X*A .+ B
    
     # Actual value of Q
     Q .= X*XAB .+ C

     # Norm of Q
     α = norm(Q)
      
     # Stop criteria
     if α<=tol
         verbose && println("Tolerance $α achieved in  $iter iterations")
         cflag = true
         break
     end
    
     # Solver the Sylvester equation
     E .= gsylv(X,A,1.0,XAB,-Q)

     # Update solution X
     X .= X .+ t*E
    
   end

   # Test for convergence
   if !cflag && verbose
      println("Warning:: Newton-Raphson iterations did not converge in $niter iterations")
   end

   # Return the solution
   return X

end



