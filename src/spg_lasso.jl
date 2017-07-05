export spg_lasso

"""
SPG_LASSO  Solve the LASSO problem\n

   SPG_LASSO is designed to solve the LASSO problem\n

   (LASSO)  minimize  ||AX - B||_2  subject to  ||X||_1 <= tau,\n

   where A is an M-by-N matrix, B is an M-vector, and TAU is a
   nonnegative scalar.  In all cases below, A can be an explicit M-by-N
   matrix or matrix-like object for which the operations  A*x  and  A'*y
   are defined (i.e., matrix-vector multiplication with A and its
   adjoint.)\n

   Also, A can be a function handle that points to a function with the
   signature\n

   v = A(w,mode)   which returns:\n
                                  v = A *w  if mode == 1;\n
                                  v = A'*w  if mode == 2. \n
   
   X = SPG_LASSO(A,B,TAU) solves the LASSO problem.\n

   X = SPG_LASSO(A,B,TAU,OPTIONS) specifies options that are set using
   SPGSETPARMS.\n

   [X,R,G,INFO] = SPG_LASSO(A,B,TAU,OPTIONS) additionally returns the
   residual R = B - A*X, the objective gradient G = A'*R, and an INFO
   structure.  (See SPGL1 for a description of this last output argument.)\n

   See also spgl1, spgSetParms, spg_bp, spg_bpdn.\n

   Copyright 2008, Ewout van den Berg and Michael P. Friedlander\n
   http://www.cs.ubc.ca/labs/scl/spgl1\n
   Id: spg_lasso.m 1074 2008-08-19 05:24:28Z ewout78 \n
"""
function spg_lasso(A ,b; tau::AbstractFloat=NaN ,options::spgOptions=spgOptions())
    
    spgl1(A, b, tau = tau, options = options)

end
