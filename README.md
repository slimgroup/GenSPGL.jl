# GenSPGL.jl
---
A Julia solver for large scale minimization problems. 

This code is an adaptation of Michael P. Friedlander, Ewout van den Berg, and Aleksandr Aravkin's MATLAB program [SPGL1](http://www.cs.ubc.ca/~mpf/spgl1/). 

## Installation
---
If you have a Github account, run the following from the Julia REPL:

    Pkg.clone("git@github.com:slimgroup/GenSPGL.jl.git")

Otherwise run: 

    Pkg.clone("https://github.com/slimgroup/GenSPGL.jl.git")

## Examples
---
**/src/Examples/spgdemo.jl**
* Example using both implicit and explicit **__A__**.

**/src/Examples/compare/**
* Example using a non linear forward function as **__A__**

