# Load a common workspace and solve a lasso problem
using MAT
vars = matread("../compare/vars_lasso.mat")

A = vars["A"]
b = vars["b"]
tau = vars["tau"]

