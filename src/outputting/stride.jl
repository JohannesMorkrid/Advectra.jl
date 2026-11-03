# ------------------------------------------------------------------------------------------
#                                   Stride & Storage Size                                   
# ------------------------------------------------------------------------------------------

# ----------------------------------- Sampling Strategy ------------------------------------

"""
    determine_strides(initial_samples, prob, total_storage_limit)

  Determines the per-diagnostic stride and checks that the cumulative storage requirement
  does not exceed `total_storage_limit` (if set). Prints the estimated total file size.
"""
function determine_strides(initial_samples, prob::SpectralODEProblem, total_storage_limit)
    strides = Int[]
    storage_requirements = Int[]
    total_storage_requirement = 0

    for (recipe, sample) in zip(prob.diagnostic_recipes, initial_samples)
        @unpack stride, storage_limit = recipe
        context = "($(recipe.method)) "

        N_samples, stride = determine_sampling_strategy(sample, stride, storage_limit, prob;
            context=context)
        storage_requirement = compute_storage_need(N_samples, stride, sample; context)
        # Accumulate
        total_storage_requirement += storage_requirement
        push!(strides, stride)
        push!(storage_requirements, storage_requirement)
    end

    # Nice to have
    println("Estimated filesize: ", format_bytes(total_storage_requirement))

    if !isempty(total_storage_limit) &&
        parse_storage_limit(total_storage_limit) < total_storage_requirement
        # TODO add nice error message showing what each diagnostic requires
        error("The Output requires $(format_bytes(total_storage_requirement)), which is \
               more than the storage_limit: $total_storage_limit.")
    end

    return strides
end

"""
    determine_sampling_strategy(sample, stride, storage_limit, prob; context="")

  Returns `(N_samples, stride)` for a single diagnostic. Interprets `stride == -1` as
  "let the program decide" based on `storage_limit`. If neither is set, every step is stored.
"""
function determine_sampling_strategy(sample, stride::Int, storage_limit,
    prob::SpectralODEProblem; context="")
    N_steps = compute_number_of_steps(prob)

    # Handle all the different scenarios, stride = -1 => let the program decide
    if !isempty(storage_limit)
        storage_bytes = parse_storage_limit(storage_limit)
        if stride == -1
            stride = recommend_stride(storage_bytes, N_steps, sample; context=context)
        else
            check_storage_size(storage_bytes, N_steps, stride, sample; context=context)
        end
    elseif stride == -1
        stride = 1
    end

    stride = validate_stride(N_steps, stride; context=context)
    N_samples = cld(N_steps, stride) + 1

    return N_samples, stride
end

"""
    compute_number_of_steps(prob)

  Number of time steps needed to integrate `prob`. The last step may be fractional and is
  rounded up because `prob.dt` is fixed.
"""
function compute_number_of_steps(prob::SpectralODEProblem)
    ceil(Int, (last(prob.tspan) - first(prob.tspan)) / prob.dt)
end

# ----------------------------------- Stride Validation ------------------------------------

"""
    validate_stride(N_steps, stride; context="")

  Returns a valid stride, adjusting and warning as needed:
  - Clamps `stride` to `N_steps` if it is too large.
  - Throws if `stride < 1`.
  - Warns when `stride` does not evenly divide `N_steps`.
"""
function validate_stride(N_steps::Int, stride::Int; context="")
    if N_steps < stride
        @warn context * "stride ($stride) exceeds total steps ($N_steps). \
                 Adjusting to stride = N_steps ($N_steps)."
        stride = N_steps
    end

    stride < 1 && throw(ArgumentError(context * "stride must be ≥ 1, got $stride"))

    if N_steps % stride != 0
        suggestion = nearest_divisor(N_steps, stride)
        @warn context * "stride ($stride) does not evenly divide N_steps ($N_steps). The \
          final output interval will be shorter. Consider using stride = $suggestion instead."
    end

    return stride
end

# ---------------------------------- Storage Size Checks -----------------------------------

"""
    recommend_stride(storage_limit, N_steps, sample; context="")

  Returns the smallest divisor of `N_steps` that keeps the total storage within
  `storage_limit`. Throws if imposed limit is too strict and informs about the minimum limit.
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

"""
    check_storage_size(storage_limit, N_steps, stride, sample; context="")

  Throws if the storage required for `(N_steps, stride, sample)` exceeds `storage_limit`.
"""
function check_storage_size(storage_limit::Int, N_steps::Int, stride::Int, sample;
    context="")
    min_stride = recommend_stride(storage_limit, N_steps, sample; context=context)
    storage_need = compute_storage_need(N_steps, stride, sample; context=context)

    if storage_need > storage_limit
        throw(ArgumentError(context * "The total output requires \
                              $(format_bytes(storage_need)), which exceeds the storage \
                              limit of $(format_bytes(storage_limit)). Consider increasing \
                              the `stride` (minimum recommended: $min_stride) or the `storage_limit`."))
    end
end

"""
    compute_storage_need(N_steps, stride, sample; context="")

  Bytes required to store `cld(N_steps, stride) + 1` copies of `sample`.
"""
function compute_storage_need(N_steps::Int, stride::Int, sample::AbstractArray; context="")
    stride < 1 && throw(ArgumentError(context * "stride must be ≥ 1, got $stride"))
    (cld(N_steps, stride) + 1) * length(sample) * sizeof(eltype(sample))
end

function compute_storage_need(N_steps, stride, sample::Number; context="")
    stride < 1 && throw(ArgumentError(context * "stride must be ≥ 1, got $stride"))
    (cld(N_steps, stride) + 1) * sizeof(eltype(sample))
end

function compute_storage_need(N_steps::Int, stride::Int, sample::Tuple; context="")
    stride < 1 && throw(ArgumentError(context * "stride must be ≥ 1, got $stride"))
    (cld(N_steps, stride) + 1) * length(sample) * sizeof(eltype(sample))
end

# Edge case
compute_storage_need(N_steps, stride, sample::Nothing; context="") = 0

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
