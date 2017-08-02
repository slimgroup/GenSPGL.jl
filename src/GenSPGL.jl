"""
Julia port of MPF's SPGL1
"""
module GenSPGL

# Modules
using MAT
using JOLI

# Types
include("spgOptions.jl")
include("spgExitCondition.jl")
include("spgInit.jl")
include("spgInfo.jl")

# Methods
include("spgl1.jl")
include("spglcore.jl")
include("spg_lasso.jl")
include("NormL1_project.jl")
include("NormL1_dual.jl")
include("NormL1_primal.jl")
include("oneprojector.jl")
include("spgline.jl")
include("spglinecurvy.jl")
include("NLfunForward.jl")
include("TraceNorm_project.jl")
include("TraceNorm_dual.jl")
include("TraceNorm_primal.jl")
include("afun.jl")

# Penalties
include("Penalty/funLS.jl")

# Examples
#include("Examples/spgdemo.jl")
include("Examples/compare/jl_cs.jl")
include("Examples/compare/jl_complex.jl")


end # module
