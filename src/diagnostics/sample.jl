# ------------------------------------------------------------------------------------------
#                                    Sample Diagnostics                                     
# ------------------------------------------------------------------------------------------

# ---------------------------------------- Density -----------------------------------------

sample_density(state, prob, time; kwargs...) = selectdim(state, ndims(prob.domain) + 1, 1)

function build_diagnostic(::Val{:sample_density}; kwargs...)
    Diagnostic(; name="Density",
               method=sample_density,
               metadata="Sampled density field")
end

# --------------------------------------- Vorticity ----------------------------------------

sample_vorticity(state, prob, time; kwargs...) = selectdim(state, ndims(prob.domain) + 1, 2)

function build_diagnostic(::Val{:sample_vorticity}; kwargs...)
    Diagnostic(; name="Vorticity",
               method=sample_vorticity,
               metadata="Sampled vorticity field")
end

# -------------------------------------- Temperature ---------------------------------------

function sample_temperature(state, prob, time; kwargs...)
    selectdim(state, ndims(prob.domain) + 1, 3)
end

function build_diagnostic(::Val{:sample_temperature}; kwargs...)
    Diagnostic(; name="Temperature",
               method=sample_temperature,
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

function build_diagnostic(::Val{:sample_potential}; kwargs...)
    Diagnostic(; name="Potential",
               method=sample_potential,
               metadata="Sampled potential field",
               assumes_spectral_state=true)
end