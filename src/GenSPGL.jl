"""
Julia port of MPF's SPGL1
"""
module GenSPGL

# Types
include("spglOptions.jl")

# Methods
include("spgl1.jl")

# Examples
include("/home/slim/klensink/.julia/v0.5/GenSPGL/examples/spgdemo.jl")

end # module
