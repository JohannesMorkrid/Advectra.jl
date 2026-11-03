# ------------------------------------------------------------------------------------------
#                                       IO Utilities                                        
# ------------------------------------------------------------------------------------------

# ---------------------------------- Storage Formatting ------------------------------------

"""
    parse_storage_limit(limit::Integer)
    parse_storage_limit(limit::AbstractString)

  Parses a storage limit as bytes. Accepts an integer (raw bytes) or a string with optional
  SI/binary units, e.g. `"10 MB"`, `"5.2 GiB"`, `"1024"`.
"""
parse_storage_limit(limit::Integer) = limit

function parse_storage_limit(limit::AbstractString)
    # Units dictionary
    units = Dict("B" => 1, "KB" => 10^3, "KIB" => 1024, "MB" => 10^6, "MIB" => 1024^2,
        "GB" => 10^9, "GIB" => 1024^3, "TB" => 10^12, "TIB" => 1024^4,
        "PB" => 10^15, "PIB" => 1024^5, "EB" => Int128(10^18), "EIB" => Int128(1024^6))

    # Match numeric value and optional unit
    # \s* matches 0 or more white spaces (\s) at the beginning ^ and end $
    # \d matches digits + -> 1 or more, * 0 or more
    # (?:\.\d) non-capturing group starting with a dot \., ? means matches zero or one
    # [a-zA-Z] matches all combinations of roman letters
    m = match(r"^\s*(\d*(?:\.\d*)?)\s*([a-zA-Z]*)\s*$", limit)
    m === nothing && throw(ArgumentError("Invalid storage limit format: $limit"))

    value = parse(Float64, m.captures[1])
    unit = uppercase(m.captures[2])
    unit = unit == "" ? "B" : unit # Converts no unit to mean byte

    haskey(units, unit) || throw(ArgumentError("Unknown storage unit: $unit"))

    multiplier = units[unit]
    return round(typeof(multiplier), value * multiplier)
end

"""
    format_bytes(nbytes)

  Formats a byte count as a human-readable string, e.g. `1032` → `"1.03 KB"`.
"""
function format_bytes(nbytes::Integer)
    units = ["B", "KB", "MB", "GB", "TB", "PB", "EB"]
    i, f = 1, float(nbytes)
    while f >= 1000 && i < length(units)
        f /= 1000
        i += 1
    end
    return string(round(f; digits=2), " ", units[i])
end

# ----------------------------------- Parameter Strings ------------------------------------

"""
    parameter_string(parameters)

  Formats parameters as a human-readable key=value string, e.g. `"κ=0.1, σ=0.02"`.
  Accepts `AbstractDict` or `NamedTuple`.
    parameter_string(parameters::AbstractDict)
    parameter_string(parameters::AbstractDict)
"""
function parameter_string(parameters::P) where {P<:AbstractDict}
    join([string(k, "=", v) for (k, v) in sort(collect(parameters))], ", ")
end

function parameter_string(parameters::P) where {P<:NamedTuple}
    join([string(k, "=", v) for (k, v) in sort(collect(pairs(parameters)))], ", ")
end

# ------------------------------------- NaN Detection --------------------------------------

"""
    assert_no_nan(state, time)

  Throws if `state[end]` is NaN, indicating a numerical breakdown.
"""
assert_no_nan(state::AbstractArray, time) =
    isnan(state[end]) && error("Breakdown occured at t=$time")

# TODO move to AdvectraGPUArraysExt.jl
function assert_no_nan(state::Union{GPUField,GPUState,AbstractGPUArray}, time)
    @allowscalar isnan(state[end]) && error("Breakdown occured at t=$time")
end
