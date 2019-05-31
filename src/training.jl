#Hardware sim functions
#these HRS and LRS values are a "little bit" of a guess but they fit the data well for BRSN002-003
LRS = 1e-3 #(siemens)
HRS = 1e-5
Tgate_ratio = 2
nrows = 7
ncols = 8

ncells = nrows*ncols
max_G = ncells*LRS
min_G = ncells*HRS

"""
Estimate the conductivity of a CMOR circuit given the CA state and the conductivity of its matching RRAM bank
"""
function conductivity(caState::Array{<:Real,2}, programMat::Array{<:Real,2})
    #do a more realistic estimate
    #currently the TGATE can only cut off ~1/2 the resistance matrix's conductivity
    caState = caState[2:end, :]
    modulation = zeros(7,8) .+ (1/Tgate_ratio) .+ ((1-1/Tgate_ratio) .* caState)

    conductances = programMat .* LRS .+ HRS

    return sum(modulation .* conductances)
end

function conductivity(caState::BitArray{2}, programMat::Union{BitArray{2}, Array{<:Real,2}})
    caState = caState .* 1.0

    if typeof(programMat) == BitArray{2}
        return conductivity(caState, programMat .* 1.0)
    else
        return conductivity(caState, programMat)
    end
end

"""
Return the conductivities for a number of different (decimal) input states.
"""
function conductivities(inputs::Array{Int,1}, rram::Array{<:Real,2}, rule::Int=60)
    n_inputs = length(inputs)

    #what's the conductivity for each input
    Gvals = zeros(Float64, n_inputs)
    for i in 1:n_inputs
        state = CA(rule, inputs[i], nrows)
        Gvals[i] = conductivity(state, rram)
    end

    return classes
end

"""
Return the log-separation of the classes. Smaller values mean a larger 'gap' can be created.
This value should be minimized.
"""
function classify(inputs::Array{Int,1}, classes::Array{Int,1}, rram::Array{<:Real,2}, rule::Int=60)
    if Set([1,-1]) != Set(unique(classes))
        error("Input classes must currently only belong to [-1,1].")
    end

    if size(rram) != (nrows, ncols)
        error(string("Inproperly shaped resistance matrix. Needs to be ", nrows, "x", ncols))
    end

    Gvals = conductivities(inputs, rram, rule)

    #what's the maximum conductivity of the low class and the minimum conductivity of the high class
    classLo = findall(f->f==-1, classes)
    classHi = findall(f->f==1, classes)

    supLoG = maximum(Gvals[classLo])
    infHiG = minimum(Gvals[classHi])

    #is there a positive separation between them
    #restrict to finite domain so that the optimizer doesn't get upset
    separation = max(infHiG - supLoG, 1e-6)
    norm_sep = separation/max_G

    #return the negative log of the positive separation as fitness value
    return -1*log(norm_sep)
end

"""
Attempt to find a RRAM bank state which creates a separation between the two desired clases of input.
"""
function optimize_reservoir(inputs::Array{Int,1}, classes::Array{Int,1}, rule::Int=60)
    lower = zeros(56)
    upper = ones(56)
    run = optimize(x->classify(inputs, classes, reshape(x,nrows,ncols), rule), lower, upper, rand(56))
    println(run)
    if !Optim.converged(run)
        return false
    elseif Optim.minimum(run) > 10
        println("\nNo solution.")
        return false
    else
        return reshape(Optim.minimizer(run),nrows,ncols)
    end
end
