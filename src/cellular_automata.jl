# # # # CA code # # # #
"""
Convert a Wolfram code (int < [0,255]) to a rule to generate the CA.
"""
function ruleGen(rule::Int)
    rule = digits(rule, base=2, pad=8)
    rules = Dict{Array{Bool,1}, Bool}()

    for i in 0:7
        x = digits(i, base=2, pad=3)
        #println(rule[i+1], " -> ", x)
        rules[[x[3], x[2], x[1]]] = rule[i+1]
    end

    return rules
end

"""
Given a elementary CA rule [0-255], calculate the iterations of that CA given a binary input.
Rows correspond to iterations, columns to the span of the 1-D CA.
"""
function CA(rule::Int, seed::Array{Bool,1}, iterations::Int)
    if !(rule >= 0 && rule <= 255)
        error("Incorrect rule [0-255]")
    end

    if (iterations < 1)
        error("Must have positive iterations.")
    end

    rules = ruleGen(rule)

    n = length(seed)
    result = falses(iterations+1, n)
    state = falses(3)

    result[1,:] = seed

    for it in 2:iterations+1
        for i in 1:n
            left = mod1(i-1, n)
            right = mod1(i+1, n)

            state[1] = result[it-1, left]
            state[2] = result[it-1, i]
            state[3] = result[it-1, right]

            result[it, i] = rules[state]
        end
    end

    return result
end

function CA(rule::Int, seed::BitArray{1}, iterations::Int)
    seed = convert(Array{Bool,1}, seed)

    return CA(rule, seed, iterations)
end

function CA(rule::Int, seed::Array{<:Int, 1}, iterations::Int)
    int_to_bool = x->convert(Array{Bool,1},x)
    return CA(rule, int_to_bool(seed), iterations)
end

function CA(rule::Int, seed::BitArray{1}, iterations::Int)
    seed = convert(Array{Bool,1}, seed)

    return CA(rule, seed, iterations)
end

function CA(rule::Int, seed::Int, iterations::Int)
    int_to_digits = x->digits(x, base=2, pad=8)
    return CA(rule, int_to_digits(seed), iterations)
end
