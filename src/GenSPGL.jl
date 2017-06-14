"""
Julia port of MPF's SPGL1
"""
module GenSPGL

# Types
include("spgOptions.jl")
include("spgExitCondition.jl")

# Methods
include("spgl1.jl")
include("spg_lasso.jl")

# Penalties
#DEVNOTE# Obviously generalize these paths
include("Penalty/funLS.jl")

# Examples
using Gallium
include("Examples/spgdemo.jl")

end # module
