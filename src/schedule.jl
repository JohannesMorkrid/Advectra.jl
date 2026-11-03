# Inspired by https://github.com/CliMA/Oceananigans.jl/blob/master/src/Utils/schedules.jl

abstract type AbstractSchedule end

# -------------------------------- Schedule Every Nth Step ---------------------------------

mutable struct Every <: AbstractSchedule
    stride::Int
    occurrence::Int
end

function Every(stride::Int)
    stride > 0 || throw(ArgumentError("Every: stride must be positive, got $stride"))
    Every(stride, stride)
end

isdue(schedule::Every, step, t) = step >= schedule.occurrence
advance_schedule!(schedule::Every) = schedule.occurrence += schedule.stride

function setup!(schedule::Every, Δt, t0, t; Δt_prev=nothing, kwargs...)
    if !isnothing(Δt_prev) && Δt != Δt_prev
        @warn """Every: Δt changed from $Δt_prev to $Δt on resume. \
                 Every($(schedule.stride)) will now fire every \
                 $(schedule.stride * Δt) time units instead of \
                 $(schedule.stride * Δt_prev)."""
    end
    steps = round(Int, (t - t0) / Δt)
    schedule.occurrence = (fld(steps, schedule.stride) + 1) * schedule.stride
end

# ----------------------------------- Schedule Every dt ------------------------------------

_snap_dt(dt, Δt) = round(round(dt / Δt, RoundNearestTiesUp, digits=10), RoundNearestTiesUp) * Δt

mutable struct EveryT{T} <: AbstractSchedule
    dt::T
    occurrence::T
    n::Int
    t0::T
end

function EveryT(dt::T) where T
    dt != 0 || throw(ArgumentError("EveryT: dt must be non-zero, got $dt"))
    EveryT{T}(dt, dt, 1, zero(T))
end

function EveryT(dt::T, Δt::T) where T
    dt != 0 || throw(ArgumentError("EveryT: dt must be non-zero, got $dt"))
    sign(dt) == sign(Δt) || throw(ArgumentError(
        "EveryT: dt and Δt must have the same sign, got dt=$dt, Δt=$Δt"
    ))
    snapped_dt = _snap_dt(dt, Δt)
    EveryT{T}(snapped_dt, snapped_dt, 1, zero(T))
end

isdue(schedule::EveryT, step, t) = sign(schedule.dt) * (t - schedule.occurrence) >= 0
function advance_schedule!(schedule::EveryT)
    schedule.n += 1
    schedule.occurrence = schedule.t0 + schedule.n * schedule.dt
end

function setup!(schedule::EveryT{T}, Δt, t0, t; kwargs...) where T
    Δt != 0 || throw(ArgumentError("EveryT: Δt must be non-zero, got $Δt"))
    sign(schedule.dt) == sign(Δt) || throw(ArgumentError(
        "EveryT: dt=$(schedule.dt) and Δt=$Δt must have the same sign"
    ))
    snapped_dt = _snap_dt(schedule.dt, Δt)
    if !(snapped_dt ≈ schedule.dt)
        @warn "EveryT: requested dt=$(schedule.dt) snapped to $snapped_dt to align with Δt=$Δt"
    end
    schedule.dt = snapped_dt
    schedule.t0 = t0
    schedule.n = fld(t - t0, snapped_dt) + 1
    schedule.occurrence = t0 + schedule.n * snapped_dt
end

# --------------------------- Schedule At Specific Occurrences -----------------------------

mutable struct At{T} <: AbstractSchedule
    occurrences::Vector{T}
    idx::Int
    direction::Int
end

At(t::T) where T = At([t])

function At(times::AbstractArray{T}) where T
    isempty(times) && throw(ArgumentError("At: times must not be empty"))
    is_strictly_increasing = issorted(times, lt=(<=))
    is_strictly_decreasing = issorted(times, lt=(>=))
    is_strictly_increasing || is_strictly_decreasing || throw(ArgumentError(
        "At: times must be strictly monotone (no duplicates), got $times"
    ))
    At(times, 1, is_strictly_increasing ? 1 : -1)
end

function isdue(schedule::At, step, t)
    1 <= schedule.idx <= length(schedule.occurrences) &&
        schedule.direction * (t - schedule.occurrences[schedule.idx]) >= 0
end

advance_schedule!(schedule::At) =
    1 <= schedule.idx <= length(schedule.occurrences) && (schedule.idx += 1)

function setup!(schedule::At{T}, Δt, t0, t; kwargs...) where T
    Δt != 0 || throw(ArgumentError("At: Δt must be non-zero, got $Δt"))
    if schedule.direction == 1
        sign(Δt) == 1 || throw(ArgumentError(
            "At: times must be strictly increasing for forward integration (Δt > 0), got $(schedule.occurrences)"
        ))
        endpoint = maximum(schedule.occurrences) + Δt
    else
        sign(Δt) == -1 || throw(ArgumentError(
            "At: times must be strictly decreasing for backward integration (Δt < 0), got $(schedule.occurrences)"
        ))
        endpoint = minimum(schedule.occurrences) + Δt
    end
    grid = t0:Δt:endpoint
    skipped = filter(t -> !any(g -> isapprox(t, g, atol=abs(Δt)/2), grid), schedule.occurrences)
    if !isempty(skipped)
        intervals = diff([t0; schedule.occurrences])
        suggested_Δt = reduce(
            (a, b) -> isapprox(mod(abs(a), abs(b)), 0, atol=1e-10) ? b : a,
            sort(unique(abs.(intervals)))
        )
        @warn """At: $(length(skipped)) occurrence(s) will be skipped with Δt=$Δt: \
                 $skipped. Consider Δt≈$(schedule.direction * suggested_Δt)."""
    end
    schedule.idx = schedule.direction == 1 ?
                   searchsortedlast(schedule.occurrences, t) + 1 :
                   searchsortedlast(schedule.occurrences, t, rev=true) + 1
end

# ----------------------- Schedule Builder For Use With Diagnostic -------------------------

# Constructor linker
const SCHEDULE_KWARGS = Dict(
    :every => Every,
    :dt => EveryT,
    :at => At,
)

function build_schedule(; kwargs...)
    found = filter(kv -> kv.first ∈ keys(SCHEDULE_KWARGS), kwargs)
    if length(found) > 1
        throw(ArgumentError("only one of $(join(keys(SCHEDULE_KWARGS), ", ")) can be \
              specified, got $(join(keys(found), ", "))"))
    end
    isempty(found) && return nothing
    k, v = first(found)
    return SCHEDULE_KWARGS[k](v)
end