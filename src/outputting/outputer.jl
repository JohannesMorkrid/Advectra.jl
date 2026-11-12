# ------------------------------------------------------------------------------------------
#                                         Outputer                                          
# ------------------------------------------------------------------------------------------

include("io_utilities.jl")        # parse_storage_limit, format_bytes, parameter_string, assert_no_nan
include("stride.jl")           # determine_strides, determine_sampling_strategy, divisor helpers
include("abstract_writer.jl")  # AbstractWriter + default no-op interface
include("hdf5_writer.jl")      # HDF5Writer  <: AbstractWriter
include("local_writer.jl")     # LocalWriter <: AbstractWriter  (TODO: partially implemented)
include("sampling.jl")      # initialize_diagnostics, handle_output!, sample_diagnostic!, ...
include("output.jl")           # Output struct + constructor (depends on all of the above)