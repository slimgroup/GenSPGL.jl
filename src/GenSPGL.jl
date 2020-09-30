"""
Julia port of MPF's SPGL1
"""
module GenSPGL

# Modules
using Printf
using Nullables
using LinearAlgebra
using MAT
using JLD
using Arpack
using Requires

const InTypeF = Union{AbstractArray, Function}
const InType = Union{AbstractArray}

function __init__()
    @require JOLI = "bb331ad6-a1cf-11e9-23da-9bcb53c69f6f"
    @eval using JOLI
    global const InTypeF = Union{joAbstractLinearOperator, AbstractArray, Function}
    global const InType = Union{joAbstractLinearOperator, AbstractArray}
end

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

# Composite
include("Composite/funCompR1.jl")
include("Composite/funCompR2.jl")

# Examples
#include("Examples/spgdemo.jl")
include("Examples/compare/jl_cs.jl")
include("Examples/pres/ex_cs.jl")
include("Examples/scaletest/scale_cs.jl")
include("Examples/compare/jl_complex.jl")

end # module
