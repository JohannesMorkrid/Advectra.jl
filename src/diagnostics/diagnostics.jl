# ------------------------------------------------------------------------------------------
#                                        Diagnostics                                        
# ------------------------------------------------------------------------------------------

# ---------------------------------- Diagnostic Struct -------------------------------------

"""
    Diagnostic{M, S, A, K}

  Runtime diagnostic object. Owns its schedule so `should_sample`, `advance_schedule!`
  and `setup!` dispatch directly on it.
"""
struct Diagnostic{M<:Function,S<:AbstractSchedule,A<:Tuple,K<:NamedTuple}
    name::String
    method::M
    schedule::S
    metadata::String
    assumes_spectral_state::Bool
    stores_data::Bool
    args::A
    kwargs::K

    function Diagnostic(; name::String,
        method::M,
        schedule::S=Every(1),
        metadata::String="",
        assumes_spectral_state::Bool=false,
        stores_data::Bool=true,
        args::A=(),
        kwargs::K=NamedTuple()) where {M<:Function,S<:AbstractSchedule,A<:Tuple,K<:NamedTuple}
        new{M,S,A,K}(name, method, schedule, metadata, assumes_spectral_state, stores_data,
            args, kwargs)
    end
end

"""
    (diagnostic::Diagnostic)(state, prob, time)

  Apply the diagnostic method with its stored args and kwargs.
"""
@inline function (diagnostic::Diagnostic)(state, prob, time)
    diagnostic.method(state, prob, time, diagnostic.args...; diagnostic.kwargs...)
end

should_sample(diagnostic::Diagnostic, step, t) = isdue(diagnostic.schedule, step, t)
advance_schedule!(diagnostic::Diagnostic) = advance_schedule!(diagnostic.schedule)
setup!(dia::Diagnostic, args...; kwargs...) = setup!(dia.schedule, args...; kwargs...)

function Base.show(io::IO, m::MIME"text/plain", diagnostic::Diagnostic)
    print(io, diagnostic.name, ": (spectral=", diagnostic.assumes_spectral_state,
        ", stores_data=", diagnostic.stores_data, ")")
    length(diagnostic.args) != 0 && print(io, ", args=", diagnostic.args)
    length(diagnostic.kwargs) != 0 && print(io, ", kwargs=", diagnostic.kwargs)
end

# -------------------------------- Building Of Diagnostics ---------------------------------

function build_diagnostic(method::Function; kwargs...)
    build_diagnostic(Val(Symbol(method)); kwargs...)
end

# ---------------------------------- Diagnostics Recipe ------------------------------------

"""
    DiagnosticRecipe

  Lightweight user-facing specification for a diagnostic. Holds the method, an already-
  constructed `AbstractSchedule`, an optional storage limit, and any extra kwargs passed
  to `build_diagnostic`. The schedule is constructed here (no `prob` needed); `setup!`
  is called later in `initialize_diagnostics` once `prob` and `t0` are known.
"""
struct DiagnosticRecipe
    method::Function
    schedule::AbstractSchedule
    storage_limit::String
    kwargs::NamedTuple

    function DiagnosticRecipe(method::Function; storage_limit::String="", kwargs...)
        # Default schedule: sample every step
        schedule = build_schedule(; kwargs...)
        # If no schedule kwarg and no storage_limit, default to Every(1), 
        # if a storage_limit is provided it means the program should resolve the stride.
        isnothing(schedule) && (schedule = isempty(storage_limit) ? Every(1) : unresolved_schedule())
        remaining_kwargs = NamedTuple(filter(kv -> kv.first ∉ keys(SCHEDULE_KWARGS), kwargs))
        new(method, schedule, storage_limit, remaining_kwargs)
    end
end

# ---------------------------------- Diagnostics Macro ------------------------------------

macro diagnostics(expr)
    if expr.head == :vect  # Vector literal
        return :([$(map(parse_diagnostic_expr, expr.args)...)])
    elseif expr.head != :block
        return :([$(parse_diagnostic_expr(expr))])
    end

    # Multiple expressions in a block
    recipes = Expr[]
    for line in expr.args
        line isa LineNumberNode && continue
        isnothing(line) && continue

        # Handle tuple expressions (when commas are used)
        if line isa Expr && line.head == :tuple
            for item in line.args
                push!(recipes, parse_diagnostic_expr(item))
            end
        else
            push!(recipes, parse_diagnostic_expr(line))
        end
    end

    # Return a vector expression directly
    return Expr(:vect, recipes...)
end

function parse_diagnostic_expr(expr)
    if expr isa Symbol
        return :(DiagnosticRecipe($expr))
    elseif expr isa Expr && expr.head == :call
        method = expr.args[1]
        kwargs = []
        for arg in expr.args[2:end]
            if arg isa Expr && arg.head == :parameters
                # Unwrap the parameters block (from semicolon syntax)
                append!(kwargs, arg.args)
            else
                # Regular keyword argument
                push!(kwargs, arg)
            end
        end

        return :(DiagnosticRecipe($method; $(kwargs...)))
    elseif expr.head == :(=)
        error("Aliases, alias = method(kwargs...), is not supported for diagnostics.")
    else
        error("Invalid diagnostic syntax: $expr")
    end
end

# ---------------------------- Include Implemented Diagnostics -----------------------------

include("CFL.jl")
include("COM.jl")
include("display.jl")
include("energy_integrals.jl")
include("probe.jl")
include("profiles.jl")
include("progress.jl")
include("sample.jl")
include("spectral.jl")

# ---------------------------------- Default Diagnostics -----------------------------------

const DEFAULT_DIAGNOSTICS = @diagnostics [progress]

requires_operator(::Val{method}; kwargs...) where {method} = OperatorRecipe[]

function required_operators(diagnostic_recipes::Vector{<:DiagnosticRecipe})
    vcat([requires_operator(Val(Symbol(recipe.method)); recipe.kwargs...)
          for recipe in diagnostic_recipes]...)
end