# ------------------------------------------------------------------------------------------
#                                    Sample Diagnostics                                     
# ------------------------------------------------------------------------------------------

# ---------------------------------------- Density -----------------------------------------

sample_density(state, prob, time; kwargs...) = selectdim(state, ndims(prob.domain) + 1, 1)

function build_diagnostic(::Val{:sample_density}; schedule=Every(1), kwargs...)
    Diagnostic(; name="Density",
               method=sample_density,
               schedule=schedule,
               metadata="Sampled density field")
end

# --------------------------------------- Vorticity ----------------------------------------

sample_vorticity(state, prob, time; kwargs...) = selectdim(state, ndims(prob.domain) + 1, 2)

function build_diagnostic(::Val{:sample_vorticity}; schedule=Every(1), kwargs...)
    Diagnostic(; name="Vorticity",
               method=sample_vorticity,
               schedule=schedule,
               metadata="Sampled vorticity field")
end

# -------------------------------------- Temperature ---------------------------------------

function sample_temperature(state, prob, time; kwargs...)
    selectdim(state, ndims(prob.domain) + 1, 3)
end

function build_diagnostic(::Val{:sample_temperature}; schedule=Every(1), kwargs...)
    Diagnostic(; name="Temperature",
               method=sample_temperature,
               schedule=schedule,
               metadata="Sampled temperature field")
end

# --------------------------------------- Potential ----------------------------------------

function sample_potential(state_hat, prob, time; kwargs...)
    @unpack operators, domain = prob
    @unpack solve_phi = operators
    slices = eachslice(state_hat; dims=ndims(state_hat))
    n_hat = slices[1]
    Ω_hat = slices[2]
    ϕ = bwd(domain) * solve_phi(n_hat, Ω_hat)
    return ϕ
end

requires_operator(::Val{:sample_potential}; kwargs...) = [OperatorRecipe(:solve_phi)]

function build_diagnostic(::Val{:sample_potential}; schedule=Every(1), kwargs...)
    Diagnostic(; name="Potential",
               method=sample_potential,
               schedule=schedule,
               metadata="Sampled potential field",
               assumes_spectral_state=true)
end
