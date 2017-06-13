"""
Julia port of MPF's SPGL1
"""
module GenSPGL

# Types
include("spgOptions.jl")

# Methods
include("spgl1.jl")
include("spg_lasso.jl")

# Penalties
#DEVNOTE# Obviously generalize these paths
include("/home/slim/klensink/.julia/v0.5/GenSPGL/src/Penalty/funLS.jl")

# Examples
using Gallium
include("/home/slim/klensink/.julia/v0.5/GenSPGL/examples/spgdemo.jl")

end # module
