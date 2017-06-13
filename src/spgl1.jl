# SPGL

export spgl1

"""
This will contain info on use of spgl1
"""
function spgl1(A::AbstractArray, b::AbstractArray;
                    tau::Number = NaN, sigma::Number = NaN, options::spglOptions = spgl_setparms())

# Check Tau and Sigma
if isnan(tau) & isnan(sigma)
    tau = 0.
    sigma = 0.
    singleTau = false
elseif isnan(sigma)
    singleTau = true
else
    if isnan(tau)
        tau = 0.
    end
    singleTau = false
end



end #func

