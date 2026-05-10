# ------------------------------------------------------------------------------------------
#                                        Diagnostics                                        
# ------------------------------------------------------------------------------------------

# ---------------------------------- Main Functionality ------------------------------------

struct Diagnostic{N<:AbstractString,M<:Function,L<:AbstractString,A<:Tuple,
                  K<:NamedTuple}
    name::N
    method::M
    metadata::L
    assumes_spectral_state::Bool
    stores_data::Bool
    args::A
    kwargs::K

    function Diagnostic(; name::AbstractString,
                        method::Function,
                        metadata::String="",
                        assumes_spectral_state::Bool=false,
                        stores_data::Bool=true,
                        args::Tuple=(),
                        kwargs=NamedTuple())
        new{typeof(name),typeof(method),typeof(metadata),typeof(args),
            typeof(kwargs)}(name, method, metadata, assumes_spectral_state, stores_data,
                            args, kwargs)
    end
end

"""
    (diagnostic::Diagnostic)(state, prob, time)

  Apply diagnostic `method` with the right `args` and `kwargs` stored in the `Diagnostic`.
"""
@inline function (diagnostic::Diagnostic)(state, prob, time)
    diagnostic.method(state, prob, time, diagnostic.args...; diagnostic.kwargs...)
end

function Base.show(io::IO, m::MIME"text/plain", diagnostic::Diagnostic)
    print(io, diagnostic.name, ": (spectral=", diagnostic.assumes_spectral_state,
          ", stores_data=", diagnostic.stores_data, ")")
    length(diagnostic.args) != 0 ? print(io, ", args=", diagnostic.args) : nothing
    length(diagnostic.kwargs) != 0 ? print(io, ", kwargs=", diagnostic.kwargs) : nothing
end

# -------------------------------- Building Of Diagnostics ---------------------------------

function build_diagnostic(method::Function; kwargs...)
    build_diagnostic(Val(Symbol(method)); kwargs...)
end

# ---------------------------------- Diagnostics Recipe ------------------------------------

struct DiagnosticRecipe
    method::Function
    stride::Int
    storage_limit::String
    kwargs::NamedTuple

    function DiagnosticRecipe(method::Function; stride::Int=-1, storage_limit="", kwargs...)
        new(method, stride, storage_limit, NamedTuple(kwargs))
    end
end

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

"""

  Allows for single method()
"""
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
include("fluxes.jl")
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