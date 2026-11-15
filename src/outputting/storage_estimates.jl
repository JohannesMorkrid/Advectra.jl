# ------------------------------------------------------------------------------------------
#                                   Stride & Storage Size
# ------------------------------------------------------------------------------------------

# ---------------------------------- Schedule Resolution -----------------------------------

"""
    resolve_schedules!(diagnostics, initial_samples, prob)

  Resolves any unresolved schedules in-place. For each diagnostic whose schedule is
  `Every(-1, typemax(Int))`, infers the best `Every(n)` from the per-diagnostic
  `storage_limit`. Errors if the diagnostic has neither a schedule nor a `storage_limit`.
"""
function resolve_schedules!(diagnostics, initial_samples, prob)
    N_steps = compute_number_of_steps(prob)
    t0 = first(prob.tspan)

    for (diagnostic, sample, recipe) in zip(diagnostics, initial_samples,
                                            prob.diagnostic_recipes)
        if schedule_unresolved(diagnostic.schedule)
            context = "($(diagnostic.name)) "
            storage_bytes = parse_storage_limit(recipe.storage_limit)
            n = recommend_stride(storage_bytes, N_steps, sample; context)

            @info context * "inferred Every($n) from storage_limit=$(recipe.storage_limit)"
            diagnostic.schedule.stride = n
            diagnostic.schedule.occurrence = n
        end

        # Align schedule to t0, separate call to `setup_schedules!` is need to sync to checkpoint
        setup!(diagnostic, prob.dt, t0, t0)
    end
end

# --------------------------------- Storage Estimation ------------------------------------

"""
    StorageEstimate

  Per-diagnostic storage estimate produced by `storage_estimates`.
"""
struct StorageEstimate
    name::String
    schedule::String
    sample_size::String
    n_samples::Int
    bytes::Int
end

"""
    storage_estimates(diagnostics, initial_samples, prob)

  Estimates storage requirements for each diagnostic based on its (now resolved) schedule
  and sample size. Pure computation, no side effects.
"""
function storage_estimates(diagnostics, initial_samples, prob)
    N_steps = compute_number_of_steps(prob)
    estimates = StorageEstimate[]

    for (diagnostic, sample) in zip(diagnostics, initial_samples)
        N_samples = estimated_samples(diagnostic.schedule, N_steps, prob)
        bytes = N_samples * sample_bytes(sample)
        push!(estimates,
              StorageEstimate(diagnostic.name,
                              string(diagnostic.schedule),
                              sample_size_str(sample),
                              N_samples,
                              bytes))
    end

    return estimates
end

# --------------------------------- Storage Reporting -------------------------------------

"""
    print_storage_report(estimates, writers)

  Prints a formatted storage estimate table to stdout.
"""
function print_storage_report(estimates::Vector{StorageEstimate}, writers)
    col_w = (20, 16, 20, 10, 14)
    sep = "─" ^ (sum(col_w) + length(col_w) * 2)
    total = sum(e.bytes for e in estimates)

    println("\nEstimated storage requirements:")
    println("  ", rpad("Diagnostic", col_w[1]),
            rpad("Schedule", col_w[2]),
            rpad("Sample size", col_w[3]),
            rpad("N samples", col_w[4]),
            rpad("Uncompressed", col_w[5]))
    println("  ", sep)
    for e in estimates
        println("  ", rpad(e.name, col_w[1]),
                rpad(e.schedule, col_w[2]),
                rpad(e.sample_size, col_w[3]),
                rpad(string(e.n_samples), col_w[4]),
                rpad(format_bytes(e.bytes), col_w[5]))
    end
    println("  ", sep)
    println("  Total (uncompressed): ", format_bytes(total))

    # Available storage from first writer (disk space for HDF5, free RAM for LocalWriter)
    if !isempty(writers)
        avail = available_storage(first(writers))
        if !isnothing(avail)
            println("  Available storage:  ", format_bytes(avail),
                    "  (note: HPC quota limits may apply — set storage_limit explicitly)")
        end
    end
    println("  Compression note: blosc typically reduces field data by 2-5x\n")
end # Written by Claude

# --------------------------------- Storage Validation ------------------------------------

"""
    validate_storage!(estimates, prob, writers, total_storage_limit)

  Checks per-diagnostic and total storage estimates against limits. Errors if any limit
  is exceeded. Uses available disk space as the default total limit when `total_storage_limit`
  is empty. Pass `storage_limit=nothing` on `Output` to skip the disk space check.
"""
function validate_storage!(estimates::Vector{StorageEstimate}, prob, writers,
                           total_storage_limit)

    # Per-diagnostic limits
    for (estimate, recipe) in zip(estimates, prob.diagnostic_recipes)
        isempty(recipe.storage_limit) && continue
        limit_bytes = parse_storage_limit(recipe.storage_limit)
        estimate.bytes > limit_bytes &&
            error("Diagnostic \"$(estimate.name)\" requires $(format_bytes(estimate.bytes)), " *
                  "exceeding its storage_limit=$(recipe.storage_limit).")
    end

    # Explicit, disk/RAM available, or skip
    has_explicit_limit = !isnothing(total_storage_limit) && !isempty(total_storage_limit)
    limit_bytes = has_explicit_limit ? parse_storage_limit(total_storage_limit) :
                  !isempty(writers) ? available_storage(first(writers)) : nothing
    isnothing(limit_bytes) && return

    total = sum(e.bytes for e in estimates)
    total > limit_bytes &&
        error("Total estimated storage $(format_bytes(total)) exceeds " *
              (has_explicit_limit ? "storage_limit=$total_storage_limit" :
               "available storage") *
              " ($(format_bytes(limit_bytes))). Reduce sampling frequency or " *
              "set storage_limit=nothing to skip this check.")
end

# ----------------------------------- Schedule Helpers -------------------------------------

"""
    estimated_samples(schedule, N_steps, prob)

  Estimates the number of samples a schedule will produce over `N_steps`.
"""
estimated_samples(schedule::Every, N_steps, prob) = cld(N_steps, schedule.stride) + 1
function estimated_samples(schedule::EveryT, N_steps, prob)
    floor(Int, (last(prob.tspan) - first(prob.tspan)) / schedule.dt) + 1
end
estimated_samples(schedule::At, N_steps, prob) = length(schedule.occurrences)

"""
    sample_bytes(sample)

  Uncompressed bytes for a single sample.
"""
sample_bytes(sample::AbstractArray) = length(sample) * sizeof(eltype(sample))
sample_bytes(sample::Number) = sizeof(eltype(sample))
sample_bytes(sample::Tuple) = length(sample) * sizeof(eltype(sample))
sample_bytes(::Nothing) = 0

"""
    sample_size_str(sample)

  Human-readable sample dimensions + element type, e.g. `"512×512 Float32"`.
"""
function sample_size_str(sample::AbstractArray)
    join(size(sample), "×") * " " * string(eltype(sample))
end
sample_size_str(sample::Number) = "scalar " * string(typeof(sample))
sample_size_str(sample::Tuple) = "$(length(sample))-tuple " * string(eltype(sample))
sample_size_str(::Nothing) = "N/A"

"""
    available_disk_space(filepath)

  Returns available bytes on the filesystem containing `filepath`, or `nothing` if the
  query fails (e.g. on some network filesystems).
"""
function available_disk_space(filepath::AbstractString)
    try
        return Base.Filesystem.diskstat(dirname(abspath(filepath))).available
    catch
        return nothing
    end
end

"""
    compute_number_of_steps(prob)

  Number of time steps needed to integrate `prob`. The last step may be fractional and is
  rounded up because `prob.dt` is fixed.
"""
function compute_number_of_steps(prob::SpectralODEProblem)
    ceil(Int, (last(prob.tspan) - first(prob.tspan)) / prob.dt)
end

# ---------------------------------- Storage Size Checks -----------------------------------

"""
    recommend_stride(storage_limit, N_steps, sample; context="")

  Returns the smallest divisor of `N_steps` that keeps the total storage within
  `storage_limit`. Throws if the imposed limit is too strict.
"""
function recommend_stride(storage_limit::Int, N_steps::Int, sample::AbstractArray;
                          context="")
    field_bytes = length(sample) * sizeof(eltype(sample))
    # Determine how many fields can be fully stored
    max_samples = storage_limit ÷ field_bytes

    # At least two should be stored, start and end
    if max_samples < 2
        throw(ArgumentError(context * "The storage limit ($(format_bytes(storage_limit))) \
          is too small. The sample alone requires $(format_bytes(field_bytes)), however at \
          least two samples are required (minimum limit: $(format_bytes(2*field_bytes)))."))
    end

    # Compute minimum stride to achieve the max number of samples
    min_stride = ceil(Int, N_steps / (max_samples - 1))
    # Return closest divisor that does not exceed storage limit, while N_steps ≥ min_stride
    return next_divisor(N_steps, min_stride)
end

# Edge case
recommend_stride(storage_limit, N_steps, sample::Nothing; context="") = 1

# ------------------------------------ Divisor Helpers -------------------------------------

"""
    next_divisor(N, target)

  Smallest divisor of `N` that is ≥ `target`. Returns `N` if none exists.
"""
function next_divisor(N::Int, target::Int)
    recommended = N
    for d in 1:floor(Int, sqrt(N))
        if N % d == 0
            for c in (d, N ÷ d)
                c >= target && c < recommended && (recommended = c)
            end
        end
    end
    return recommended
end

"""
    nearest_divisor(N, target)

  Divisor of `N` closest to `target`. Returns `1` if `N` is prime.
"""
function nearest_divisor(N::Int, target::Int)
    closest, smallest = 1, abs(target - 1)
    for d in 1:floor(Int, sqrt(N))
        if N % d == 0
            for c in (d, N ÷ d)
                dist = abs(c - target)
                if dist < smallest || (dist == smallest && c < closest)
                    closest, smallest = c, dist
                end
            end
        end
    end
    return closest
end
