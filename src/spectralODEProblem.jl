
# Inspired by SciMLBase, singleton type
"""
    Nullparamters()
  Singleton to signalize that the user has not specified any parameters.
"""
struct NullParameters end

abstract type AbstractODEProblem{isinplace} end

"""
    SpectralODEProblem(L::Function, N::Function, u0, domain::AbstractDomain, tspan;
        p=NullParameters(), dt=0.01, remove_modes::Function=remove_nothing, kwargs...)
    SpectralODEProblem(N::Function, u0, domain::AbstractDomain, tspan;
        p=NullParameters(), dt=0.01, remove_modes::Function=remove_nothing, kwargs...)
  
  Collection of data needed to specify the spectral ODE problem to be solved. The user needs
   to specify the non-linear operator `N`, with the linear operator `L` being optional and 
   otherwise assumed tp be zero. In addition the `domain`, initial condition `u0` and
   timespan `tspan` needs to be specified. The parameters `p` will be passed onto the 
   RHS/operators and the timestep `dt` is used by the temporal scheme. There is also the 
   option to add a method `remove_modes` to remove certain modes after each timestep. Other
   `kwargs` can be stored in the struct, however these are currently unused.
"""
mutable struct SpectralODEProblem{LType<:Function,NType<:Function,
                                  u0Type<:AbstractArray,u0_hatType<:AbstractArray,
                                  D<:AbstractDomain,tType,pType,
                                  operatorsType<:NamedTuple,
                                  diagnosticRecipesType<:Vector,N<:Number,
                                  RMType<:Function,K,iip} <: AbstractODEProblem{iip}
    L::LType
    N::NType
    u0::u0Type
    u0_hat::u0_hatType
    domain::D
    tspan::tType
    p::pType

    operators::operatorsType
    diagnostic_recipes::diagnosticRecipesType
    dt::N
    remove_modes::RMType
    kwargs::K

    function SpectralODEProblem(NonLinear::Function, u0, domain::AbstractDomain, tspan;
                                p=NullParameters(), dt=0.01,
                                remove_modes::Function=remove_nothing, kwargs...)

        # If no linear operator given, assume there is non and match signature
        isinplace(NonLinear) isa Val{true} ? L(du, u, operators, p, t) = (du .= zero(u)) :
        L(u, operators, p, t) = zero(u)

        SpectralODEProblem(L, NonLinear, u0, domain, tspan; p=p, dt=dt,
                           remove_modes=remove_modes, kwargs...)
    end

    function SpectralODEProblem(Linear::Function, NonLinear::Function, u0,
                                domain::AbstractDomain, tspan;
                                p=NullParameters(), dt::Number=0.01,
                                operators::Symbol=:default,
                                aliases::Vector{Pair{Symbol,Symbol}}=Pair{Symbol,Symbol}[],
                                additional_operators::Vector{<:OperatorRecipe}=OperatorRecipe[],
                                remove_modes::Function=remove_nothing,
                                diagnostics::Vector{<:DiagnosticRecipe}=DiagnosticRecipe[],
                                kwargs...)

        # Prepare data structures
        u0 = prepare_initial_condition(u0, domain)
        u0_hat = prepare_spectral_coefficients(u0, domain)

        # Remove unwanted modes
        remove_modes(u0_hat, domain)

        # Handle time related things
        length(tspan) != 2 ? throw("tspan should have exactly two elements") : nothing
        tspan = convert.(domain.precision, promote(first(tspan), last(tspan)))
        dt = convert(domain.precision, dt)

        # Returns a NamedTuple with `SpectralOperator`s
        ops = build_operators(domain; operators, aliases, additional_operators, diagnostics,
                              kwargs...)

        # Makes the rhs follow the signature used by SciML
        L, N = prepare_functions(Linear, NonLinear, ops)

        new{typeof(L),typeof(N),typeof(u0),typeof(u0_hat),typeof(domain),typeof(tspan),
            typeof(p),typeof(ops),typeof(diagnostics),typeof(dt),typeof(remove_modes),
            typeof(kwargs),isinplace(Linear, NonLinear)}(L, N, u0, u0_hat, domain, tspan, p,
                                                         ops, diagnostics, dt, remove_modes,
                                                         kwargs)
    end
end

# --------------------------------- Construction Related -----------------------------------

"""
"""
function prepare_functions(Linear::Function, NonLinear::Function, operators::NamedTuple)
    if isinplace(Linear, NonLinear) isa Val{true}
        L(du, u, p, t) = Linear(du, u, operators, p, t)
        N(du, u, p, t) = NonLinear(du, u, operators, p, t)
    else
        L = (u, p, t) -> Linear(u, operators, p, t)
        N = (u, p, t) -> NonLinear(u, operators, p, t)
    end
    return L, N
end

"""
"""
function prepare_initial_condition(u0, domain::Domain)
    # Transform to the correct MemoryType for physical space
    u0 |> domain.MemoryType{eltype(fwd(domain))}
end

# -------------------------- Spectral Coefficent Initialization ----------------------------

"""
"""
function prepare_spectral_coefficients(u0, domain::Domain)
    # Allocate data structure for spectral modes
    u0_hat = allocate_coefficients(u0, domain)

    # Compute spectral initial conditions
    spectral_transform!(u0_hat, get_fwd(domain), u0)

    return u0_hat
end

"""
    allocate_coefficients(u0, transformplans::AbstractTransformPlans)
  
  Recursively iterates trough the initial condition data structure to try to get to the 
  lowest level Array and then allocates the needed shape after applying the fwd transform.
"""
allocate_coefficients(u0, domain::Domain) = _allocate_coefficients(u0, domain)

# TODO perhaps clean up this logic
function _allocate_coefficients(u0::AbstractArray{<:Number}, domain::Domain)
    # Allocate array for spectral modes 
    sz = size(get_bwd(domain))
    allocation_size = (sz..., size(u0)[(length(sz)+1):end]...)
    return zeros(eltype(get_bwd(domain)), allocation_size) |> domain.MemoryType
end

function _allocate_coefficients(u0::AbstractArray{<:AbstractArray}, domain::Domain)
    [_allocate_coefficients(u, domain) for u in u0]
end

# --------------------------------- Building Of Operators ----------------------------------

"""
    prepare_operator_recipes(operators::Symbol,
                                  additional_operators::Vector{<:OperatorRecipe},
                                  diagnostics::Vector{<:DiagnosticRecipe})

  Use switches to get a list of `OperatorRecipe`s and append the `additional_operators` and
  `required_operators` by the `diagnostics` to it.
"""
function prepare_operator_recipes(operators::Symbol,
                                  additional_operators::Vector{<:OperatorRecipe},
                                  diagnostics::Vector{<:DiagnosticRecipe})
    # Determine operators trough switches
    recipes = get_operator_recipes(operators)
    # Determine which operators are required by the diagnostics
    required = required_operators(diagnostics)
    # Combine with additional_operators 
    return vcat(recipes, additional_operators, required)
end

function ensure!(cache, recipe, domain, problem_kwargs)
    # To not construct the same operator twice
    if haskey(cache, recipe)
        return cache[recipe]
    end

    # Makes the code more readable
    operator = recipe.op

    # Get dependencies and construct recursively
    dependencies = [dependency.op => ensure!(cache, dependency, domain, problem_kwargs)
                    for dependency in operator_dependencies(Val(operator), Domain)]

    # Collect kwargs from recipe, problem and combine with dependency references
    kwargs = (; recipe.kwargs..., problem_kwargs..., dependencies...)

    # Build operator and add to cache
    cache[recipe] = build_operator(Val(operator), domain; kwargs...)
end

# TODO fix aliases method issue [#27](https://github.com/JohannesMorkrid/Advectra.jl/issues/27)
function add_aliases!(operators, aliases, cache)
    #for alias in aliases
    #    # Get last, as thats whats being aliases
    #    operator = last(alias)
    #    cache[]
    #end
    #vcat()
    #[first(alias) => cache[last(aliases)] for alias in aliases]
    return operators
end

function build_operators(domain::AbstractDomain; operators::Symbol=:default,
                         aliases::Vector{Pair{Symbol,Symbol}}=Pair{Symbol,Symbol}[],
                         additional_operators::Vector{<:OperatorRecipe}=OperatorRecipe[],
                         diagnostics::Vector{<:DiagnosticRecipe}=DiagnosticRecipe[],
                         problem_kwargs...)
    # Collects all recipes needed to be built
    recipes = prepare_operator_recipes(operators, additional_operators, diagnostics)

    # Used to only have to build each operator once
    cache = Dict{OperatorRecipe,SpectralOperator}()

    # Construct NamedTuple
    spectral_operators = (;
                          [recipe.alias => ensure!(cache, recipe, domain, problem_kwargs)
                           for recipe in recipes]...)

    # Sort out aliases
    spectral_operators = add_aliases!(spectral_operators, aliases, cache)

    return spectral_operators
end

# ---------------------------------------- Helpers -----------------------------------------

spectral_size(prob::SpectralODEProblem) = size(prob.u0_hat)

get_precision(prob::SpectralODEProblem) = prob.domain.precision
get_fwd(prob::SpectralODEProblem) = get_fwd(prob.domain)
get_bwd(prob::SpectralODEProblem) = get_bwd(prob.domain)

function get_kwargs(prob::SpectralODEProblem)
    (; [field => getfield(prob, field) for field in fieldnames(SpectralODEProblem)]...)
end

function build_diagnostic(method::Function, prob::SpectralODEProblem; kwargs...)
    prob_kwargs = get_kwargs(prob)
    build_diagnostic(Val(Symbol(method)); prob_kwargs..., kwargs...)
end

"""
    isinplace(prob::AbstractODEProblem{iip}) where {iip}
    isinplace(f::Function)
    isinplace(L::Function, N::Function)

  Checks whether or not the functions are in place. Works like a trait.
"""
isinplace(prob::AbstractODEProblem{iip}) where {iip} = iip

function isinplace(L::Function, N::Function)
    inplace = isinplace(L)
    inplace == isinplace(N) ? inplace : error("Mismatch in function signatures: Both \
    `prob.L` and `prob.N` must have the same signature. `prob.L` is 
    $(inplace ? "in-place" : "out-of-place"), while `prob.N` is 
    $(!inplace ? "in-place" : "out-of-place").")
end

# Inspired by https://github.com/SciML/SciMLBase.jl/blob/d1072adfcb061db6617972d4d5b2b6610ab32839/src/utils.jl#L6
function isinplace(f::Function)
    nargs = [m.nargs - 1 for m in methods(f)]
    inplace = any(x -> x == 5, nargs)
    outofplace = any(x -> x == 4, nargs)

    if inplace
        return Val(true)
    elseif outofplace
        return Val(false)
    end

    error("The function `$f` must have a valid signature.
    Expected either:
        - In-place: `$f(du, u, d, p, t)` (5 arguments, modifies `du` in-place), or
        - Out-of-place: `$f(u, d, p, t)` (4 arguments, returns a new value).
    However, no methods of `$f` match these signatures.")
end

get_linear_operator(prob::SpectralODEProblem) = prob.L.Linear
get_nonlinear_operator(prob::SpectralODEProblem) = prob.N.NonLinear

function Base.show(io::IO, m::MIME"text/plain", prob::SpectralODEProblem)
    print(io, nameof(typeof(prob)), "(", nameof(prob.L), ",", nameof(prob.N), ";dt=",
          prob.dt)
    typeof(prob.p) == NullParameters ? nothing : print(io, ",p=", prob.p)
    println(io, "):")
    println(io, "in-place: ", isinplace(prob) isa Val{true})
    println(io, "remove_modes: ", nameof(prob.remove_modes))
    print(io, "domain: ")
    show(io, m, prob.domain)
    print(io, "\ntimespan: ")
    show(io, prob.tspan)
    print(io, "\nu0: ")
    show(io, m, prob.u0)
end
