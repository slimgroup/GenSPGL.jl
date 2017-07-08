"""
Julia port of MPF's SPGL1
"""
module GenSPGL

# Modules
using MAT

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
include("spglinecurvy.jl")
include("NLfunForward.jl")

# Penalties
#DEVNOTE# Obviously generalize these paths
include("Penalty/funLS.jl")

# Examples
include("Examples/spgdemo.jl")
include("../compare/jl_lasso.jl")
end # module
