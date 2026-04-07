# ------------------------------------------------------------------------------------------
#                                    Display Diagnostic                                     
# ------------------------------------------------------------------------------------------

function plot_field(domain, field, time=-1; field_name="", digits=2, kwargs...)
    ending = (time == -1 ? ")" : ", t = $(round(time, digits=digits)))")
    title = field_name * "(x, y" * ending
    display(heatmap(domain, real(field); aspect_ratio=:equal, xlabel="x", ylabel="y",
                    title=title, kwargs...))
end

# ---------------------------------------- Density -----------------------------------------

function plot_density(state, prob, time; digits=2, kwargs...)
    n = selectdim(state, ndims(prob.domain) + 1, 1) |> Array
    plot_field(prob.domain, n, time; field_name=L"n", digits=digits, kwargs...)
end

function build_diagnostic(::Val{:plot_density}; dt, kwargs...)
    kwargs = (; digits=ceil(Int, -log10(dt)))
    Diagnostic(; name="Plot density",
               method=plot_density,
               metadata="Display density",
               stores_data=false,
               kwargs=kwargs)
end

# --------------------------------------- Vorticity ----------------------------------------

function plot_vorticity(state, prob, time; digits=2, kwargs...)
    Ω = selectdim(state, ndims(prob.domain) + 1, 2) |> Array
    plot_field(prob.domain, Ω, time; field_name=L"\Omega",
               digits=digits, color=:jet, kwargs...)
end

function build_diagnostic(::Val{:plot_vorticity}; dt, kwargs...)
    kwargs = (; digits=ceil(Int, -log10(dt)))
    Diagnostic(; name="Plot vorticity",
               method=plot_vorticity,
               metadata="Display vorticity",
               stores_data=false,
               kwargs=kwargs)
end

# -------------------------------------- Temperature ---------------------------------------

function plot_temperature(state, prob, time; digits=2, kwargs...)
    T = selectdim(state, ndims(prob.domain) + 1, 3) |> Array
    plot_field(prob.domain, T, time; field_name=L"T", digits=digits, color=:jet, kwargs...)
end

function build_diagnostic(::Val{:plot_temperature}; dt, kwargs...)
    kwargs = (; digits=ceil(Int, -log10(dt)))
    Diagnostic(; name="Plot temperature",
               method=plot_temperature,
               metadata="Display temperature",
               stores_data=false,
               kwargs=kwargs)
end

# --------------------------------------- Potential ----------------------------------------

function plot_potential(state_hat, prob, time; digits=2, kwargs...)
    @unpack operators, domain = prob
    @unpack solve_phi = operators
    slices = eachslice(state_hat; dims=ndims(state_hat))
    n_hat = slices[1]
    Ω_hat = slices[2]
    ϕ = bwd(domain) * solve_phi(n_hat, Ω_hat) |> Array
    plot_field(domain, ϕ, time; field_name=L"\phi", digits=digits, kwargs...)
end

requires_operator(::Val{:plot_potential}; kwargs...) = [OperatorRecipe(:solve_phi)]

function build_diagnostic(::Val{:plot_potential}; dt, kwargs...)
    kwargs = (; digits=ceil(Int, -log10(dt)))
    Diagnostic(; name="Display potential",
               method=plot_potential,
               metadata="Displays the potential",
               assumes_spectral_state=true,
               stores_data=false,
               kwargs=kwargs)
end