module CMOR

using PyPlot, Optim

include("src/cellular_automata.jl")
include("src/training.jl")
include("src/testing.jl")

export CA,
conductivity,
conductivities,
classify,
optimize_reservoir,
generate_config,
read_reservoir,
generate_read_program,
generate_forming_program,
read_voltages,
plot_results,
read_currents,
plot_results_hwonly

end
