# # # # CA code # # # #

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

# # # # CMOR Simulation # # # #

std_input = [1,0,0,0, 0,0,0,0]

"""
Calculate the total reservoir resistance given a CA state or input/rule and resistance matrix.
"""
function resistance(caState::BitArray{2}, resistanceMat::Array{<:Real,2})
    for i in resistanceMat
        if i <= 0.0
            error("Resistance matrix must be on positive domain")
        end
    end

    if size(caState) != size(resistanceMat)
        error("CA State and Resistance state have size mismatch.")
    end

    invsum = Float64(0.0)

    for i in 1:length(resistanceMat)
        #is this 1T1R enabled by the CA?
        if caState[i]
            #if so, its conductivity can contribute
            invsum += 1/resistanceMat[i]
        end
    end

    if invsum != 0
        return 1/invsum
    else
        return Inf
    end
end

function resistance(rule::Int, input::Union{Array{Union{Bool, Int},1}, BitArray{1}}, resistanceMat::Array{<:Real,2})
    caState = CA(rule, input, 7)

    return resistance(caState, resistanceMat)
end

"""
Given a 7x8 resistance matrix, check how good it is at performing XOR
"""
function test_XOR_resistances(rmat::Array{<:Real,2})
    null = falses(8); a = falses(8); b = falses(8); ab = falses(8)
    a[1] = true
    b[2] = true
    ab[1:2] .= true

    rule = 60
    ca_null = CA(60, null, 7)
    ca_a = CA(60, a, 7)
    ca_b = CA(60, b, 7)
    ca_ab = CA(60, ab, 7)

    return [resistance(ca_null, rmat), resistance(ca_a, rmat), resistance(ca_b, rmat), resistance(ca_ab, rmat)]
end


"""
Check how good a resistance matrix is at being a linear classifier (separating 3 high from 3 low values).
"""
function linear_fitness(rule::Int, r::Array{<:Real,1})
    n = 8

    if length(r) != n^2
        error("Wrong resistance matrix length.")
    end

    #set up the low values
    a = [falses(n) for i in 1:3]
    a[1][end] = true
    a[2][end-1] = true
    a[3][end-1:end] .= true

    #set up the high values
    b = [trues(n) for i in 1:3]
    b[1][end] = false
    b[2][end-1] = false
    b[3][end-1:end] .= false

    ca_a = [CA(rule, a[i], n-1) for i in 1:3]
    ca_b = [CA(rule, b[i], n-1) for i in 1:3]

    r_a = [resistance(ca_a[i], reshape(r,n,n)) for i in 1:3]
    r_b = [resistance(ca_b[i], reshape(r,n,n)) for i in 1:3]

    distance = 0.0
    for i in 1:3, j in 1:3
        distance += abs(r_a[i] - r_b[j])
    end

    return -1*distance
end

"""
Print the resistance values for the high and low examples in the linear classifier.
"""
function linear_verify(rule::Int, r::Array{<:Real,1})
    n = 8

    if length(r) != n^2
        error("Wrong resistance matrix length.")
    end

    #set up the low values
    a = [falses(n) for i in 1:3]
    a[1][end] = true
    a[2][end-1] = true
    a[3][end-1:end] .= true

    #set up the high values
    b = [trues(n) for i in 1:3]
    b[1][end] = false
    b[2][end-1] = false
    b[3][end-1:end] .= false

    ca_a = [CA(rule, a[i], n-1) for i in 1:3]
    ca_b = [CA(rule, b[i], n-1) for i in 1:3]

    r_a = [resistance(ca_a[i], reshape(r,n,n)) for i in 1:3]
    r_b = [resistance(ca_b[i], reshape(r,n,n)) for i in 1:3]

    return (r_a, r_b)
end

# # # # CMOR Testing # # # #

function check_bounds(rs::Int, cs::Int)
    if !(cs in 0:7)
        error("Column needs to be in [0,7].")
    end

    if !(rs in 1:7)
        error("Row needs to be in [1,7].")
    end
end

"""
Generates the 24-pin matrix connection string to read out the state of a CA cell given an 8-bit input.
"""
function generate_config(input::Array{<:Union{Int,Bool},1}, rs::Int, cs::Int)
    check_bounds(rs, cs)

    if length(input) != 8
        error("Need 8-bit input.")
    end

    int_to_bool = x->convert(Array{Bool,1},x)

    if typeof(input[1]) <: Int
        input = int_to_bool(input)
    end

    cs_b = int_to_bool(digits(cs, base=2, pad=3))
    rs_b = int_to_bool(digits(rs, base=2, pad=3))

    #wire GND, R-D, R-S, RVT, VDDH, unwired"
    config = "A, A, A, C, C, A, "

    #set the column selected
    for (i,x) in enumerate(cs_b)
        if x
            config *= "B, "
        else
            config *= "A, "
        end
    end

    #set the row selected
    for (i,x) in enumerate(rs_b)
        if x
            config *= "B, "
        else
            config *= "A, "
        end
    end

    #wire GND, VDDL, unused, CA-O

    config *= "A, B, A, D, "

    #set the bit input
    for (i,x) in enumerate(input)
        if x
            config *= "B, "
        else
            config *= "A, "
        end
    end

    config = config[1:end-2]

    return config
end

"""
Generates the lines to test and read out the state of the reservoir given a fixed input ([1,0,0,0, 0,0,0,0]).
"""
function read_reservoir()
    matrix = ""
    for r in 1:7, c in 0:7
        matrix = matrix * generate_config([1,0,0,0, 0,0,0,0], r, c) * "\n"
    end
    return matrix
end

"""
Generates the lines to enable a single memristor and form it.
"""
function generate_forming_program(rs::Int, cs::Int)
    check_bounds(rs, cs)

    int_to_bool = x->convert(Array{Bool,1},x)

    cs_b = int_to_bool(digits(cs, base=2, pad=3))
    rs_b = int_to_bool(digits(rs, base=2, pad=3))

    #wire GND, R-D, R-S, RVT, VDDH, unwired"
    config = "A, A, A, D, C, A, "

    #set the column selected
    for (i,x) in enumerate(cs_b)
        if x
            config *= "B, "
        else
            config *= "A, "
        end
    end

    #set the row selected
    for (i,x) in enumerate(rs_b)
        if x
            config *= "B, "
        else
            config *= "A, "
        end
    end

    #wire GND, VDDL, unused, CA-O

    config *= "A, B, A, F, "

    #set the bit input
    for i in 1:8
        config *= "A, "
    end

    config = config[1:end-2]

    return config
end

"""
Read the resistance of a single memristor, enabling its read resistor only.
"""
function generate_read_program(rs::Int, cs::Int)
    check_bounds(rs, cs)

    int_to_bool = x->convert(Array{Bool,1},x)

    cs_b = int_to_bool(digits(cs, base=2, pad=3))
    rs_b = int_to_bool(digits(rs, base=2, pad=3))

    #wire GND, R-D, R-S, RVT, VDDH, unwired"
    config = "A, A, E, C, C, A, "

    #set the column selected
    for (i,x) in enumerate(cs_b)
        if x
            config *= "B, "
        else
            config *= "A, "
        end
    end

    #set the row selected
    for (i,x) in enumerate(rs_b)
        if x
            config *= "B, "
        else
            config *= "A, "
        end
    end

    #wire GND, VDDL, unused, CA-O

    config *= "A, B, A, F, "

    #set the bit input
    for i in 1:8
        config *= "A, "
    end

    config = config[1:end-2]

    return config

end


# # # # Data Analysis # # # #

cmor_rules = [60, 90, 102, 105, 153, 165, 180, 195]

function readvoltages(file::String)
    voltages = []
    f = readlines(open(file,"r"))

    for line in f
        m = match(reg, line)
        if m != nothing
            vout = Meta.parse(m[4])
            append!(voltages, vout)
        end
    end

    return voltages
end

function plot_results(x::String, rule::Int)
    voltages = reshape(readvoltages(x),8,7), dims=2))
    vhigh = maximum(voltages)
    ax = subplot(121)
    ax[:set_title]("Software CA")
    imshow(rotr90(CA(rule, [1,0,0,0, 0,0,0,0], 7)))

    ax = subplot(122)
    ax[:set_title]("Hardware CA")
    imshow(rotr90(cat([vhigh,0,0,0, 0,0,0,0], voltages, cmap="bone")
    xlabel("CA Column")
    ylabel("CA Row")
    colorbar()
end
