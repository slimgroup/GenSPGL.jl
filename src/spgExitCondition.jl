#SPG Exit Condition Dict Constructor
import Base.print

export spgExitCondition, print

mutable struct spgExitCondition
    triggered::Nullable{Integer}
    conditions::Dict{Integer,String}
    info::Array{String,1}
end

"""
spgExitCondition(;trigger::Nullable{Integer} = Nullable{Integer}())\n

Use print(ec::spgExitCondition) to see error code and info.\n
"""
function spgExitCondition(;trigger::Nullable{Integer} = Nullable{Integer}())

    vals = ["EXIT_ROOT_FOUND",
            "EXIT_BPSOL_FOUND",
            "EXIT_LEAST_SQUARES",
            "EXIT_OPTIMAL",
            "EXIT_ITERATIONS",
            "EXIT_LINE_ERROR",
            "EXIT_SUBOPTIMAL_BP",
            "EXIT_MATVEC_LIMIT",
            "EXIT_ACTIVE_SET",
            "EXIT_AT_PARETO"]

    info = ["EXIT -- Found a root",
            "EXIT -- Found a BP solution",
            "EXIT -- Found a least-squares solution",
            "EXIT -- Optimal solution found",
            "ERROR EXIT -- Too many iterations",
            "ERROR EXIT -- Linesearch error", #DEVNOTE# Show Line
            "EXIT -- Found a suboptimal BP solution",
            "EXIT -- Maximum matrix-vector operations reached",
            "EXIT -- Found a possible active set",
            "EXIT -- Reached the pareto curve"]

    # Create conditions Dict
    conds = Dict(i => vals[i] for i=1:length(vals))

    # Construct
    exitcon = spgExitCondition(trigger, conds, info)
    
    return exitcon
end


#DEVNOTE# Not sure if this should be added to base, right now it requires GenSPGL.print()
"""
--------------------------------------------------------------------------------
GenSPGL\n
Use: print(ec::spgExitCondition)\n

Prints the triggered exit condition reached in spgl1.\n
--------------------------------------------------------------------------------
"""
function print(ec::spgExitCondition)

    if isnull(ec.triggered)
        println("Null Trigger. No exit condition reached in spgl1.")
    else
        val = ec.triggered.value
        reason = ec.conditions[val]
        info = ec.info[val]

        # Output
        println("""
        Exit Condition Number:          $(val)\n
        Exit Condition Triggered:       $(reason)\n
        Additional Information:         $(info)
        
        """)
    end

end
