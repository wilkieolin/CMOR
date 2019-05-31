function check_bounds(rs::Int, cs::Int)
    if !(cs in 0:7)
        error("Column needs to be in [0,7].")
    end

    if !(rs in 1:7)
        error("Row needs to be in [1,7].")
    end
end

"""
#SETUP REFERENCE
- A - 0.0 V - GND
- B - 1.2 V - VDDL
- C - 3.3 V - VDDH
- D - 1.2 V - RVT (memristor programming threshold)
- E - 0.3 V - RSRC (Current Sense)
- F - 0.0 V - RDRN (Memristor Ground)
- G - 0.0 A - VSENSE (Logic Read)
"""


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

function plot_results(x::String, rule::Int)
    voltages = shape_voltages(readvoltages(x))
    vhigh = maximum(voltages)

    fig = figure(figsize=(11,5))
    fig.suptitle(string("Rule ", rule))

    ax = subplot(121)
    ax.set_title("Software CA")
    im = ax.imshow(convert(Array{Int,2},(CA(rule, [1,0,0,0, 0,0,0,0], 7))),cmap="bone")
    ax.set_xlabel("Column")
    ax.set_xticks(collect(0:7), string.(collect(0:7)))
    ax.set_ylabel("Row")
    cbar = fig.colorbar(im)
    cbar.ax.set_ylabel("State")

    ax = subplot(122)
    ax.set_title("Hardware CA")
    im = ax.imshow(cat([vhigh,0,0,0, 0,0,0,0]', voltages, dims=1), cmap="bone")
    ax.set_xlabel("Column")
    ax.set_xticks(collect(0:7), string.(collect(0:7)))
    ax.set_ylabel("Row")
    cbar = fig.colorbar(im)
    cbar.ax.set_ylabel("Voltage")
end


function read_voltages(file::String)
    voltages = []
    f = readlines(open(file,"r"))
    reg = r"^DataValue ,\s+(\S+)\s+,\s+(\S+)\s+,\s+(\S+)\s+,\s+(\S+)"

    for line in f
        m = match(reg, line)
        if m != nothing
            vout = Meta.parse(m[4])
            append!(voltages, vout)
        end
    end

    return voltages
end

function shape_voltages(readings::Array{Any,1})
    if length(readings) != 56
        error("Not right number of readings (56 = 7x8).")
    end

    return transpose(reshape(readings,8,7))
end

function read_currents(file::String)
    voltages = []
    f = readlines(open(file,"r"))
    reg = r"^DataValue ,\s*(\S+)\s*,\s*(\S+)\s*,\s*(\S+)\s*,\s*(\S+)\s*,\s*(\S+)\s*,\s*(\S+)\s*"

    for line in f
        m = match(reg, line)
        if m != nothing
            vout = Meta.parse(m[5])
            append!(voltages, vout)
        end
    end

    return voltages
end

function plot_results_hwonly(x::String, rule::Int)
    voltages = shape_voltages(readvoltages(x))
    vhigh = maximum(voltages)

    fig = figure(figsize=(5,5))

    ax = subplot(111)
    ax.set_title("Hardware CA")
    im = ax.imshow(cat([vhigh,0,0,0, 0,0,0,0]', voltages, dims=1), cmap="bone")
    ax.set_xlabel("Column")
    ax.set_xticks(collect(0:7), string.(collect(0:7)))
    ax.set_ylabel("Row")
    cbar = fig.colorbar(im)
    cbar.ax.set_ylabel("Voltage")
    #cbar = fig[:colorbar](im)
    #cbar[:ax][:set_ylabel]("Voltage")
end
