
export NLfunForward


"""
"""
function NLfunForward(x, g, params)

    if isempty(g)
        f1 = params["A"]*x
        f2 = 0
    else
        f1 = params["A"]'*g
        f2 = g

    end

    return f1,f2
end
