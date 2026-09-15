using Test
using Advectra
import Advectra: get_cache, perform_step!, save_checkpoint!, restore_checkpoint,
    get_data, get_domain

# Advectra.jl#120: `Field`/`State` are meant to be usable as the underlying storage for the
# scheme `Cache`s (MSS1Cache, MSS2Cache, MSS3Cache, ...). schemes.jl itself needs no changes
# for this since it only relies on the generic AbstractArray interface (similar/copy/zero/
# broadcasting), which `Field`/`State` already implement -- these tests lock that in.

domain = Domain(16, 16; Lx=10, Ly=10)

function Linear!(du, u, operators, p, t)
    du .= p.ν .* operators.laplacian(u)
end

function NonLinear!(du, u, operators, p, t)
    du .= zero(u)
end

parameters = (ν=0.3,)
tspan = [0.0, 0.1]
dt = 1e-2

@testset "Scheme cache with Field matches plain Array" for scheme in (MSS1(), MSS2(), MSS3())
    u0_raw = initial_condition(gaussian, domain)
    u0_field = Field(copy(u0_raw), domain)

    prob_raw = SpectralODEProblem(Linear!, NonLinear!, u0_raw, domain, tspan;
        p=parameters, dt=dt, operators=:all)
    prob_field = SpectralODEProblem(Linear!, NonLinear!, u0_field, domain, tspan;
        p=parameters, dt=dt, operators=:all)

    cache_raw = get_cache(prob_raw, scheme)
    cache_field = get_cache(prob_field, scheme)

    @test cache_field.u isa Field
    @test get_domain(cache_field.u) === domain

    for (step, t) in enumerate(first(tspan):dt:(last(tspan)-dt))
        perform_step!(cache_raw, prob_raw, t)
        perform_step!(cache_field, prob_field, t)
    end

    @test get_data(cache_field.u) ≈ cache_raw.u
end

@testset "Scheme cache with State" for scheme in (MSS1(), MSS2(), MSS3())
    ic = State((:θ, :Ω), initial_condition(isolated_blob, domain), domain)
    prob = SpectralODEProblem(Linear!, NonLinear!, ic, domain, tspan;
        p=parameters, dt=dt, operators=:all)

    cache = get_cache(prob, scheme)
    @test cache.u isa State
    @test get_domain(cache.u) === domain

    for (step, t) in enumerate(first(tspan):dt:(last(tspan)-dt))
        perform_step!(cache, prob, t)
    end

    @test all(isfinite, get_data(cache.u))
end

# Advectra.jl#120 comment: checkpoint save/restore must round-trip Field/State cache
# fields, not just plain arrays. `adapt(typeof(field), raw)` alone drops the Field/State
# wrapper (it only adapts the underlying storage), so `restore_checkpoint` needs to rewrap.
@testset "Checkpoint round-trip with Field/State cache" begin
    @testset "$(typeof(scheme)), $wrapper" for scheme in (MSS1(), MSS2(), MSS3()),
        wrapper in (:Field, :State)

        ic = wrapper === :Field ? Field(initial_condition(gaussian, domain), domain) :
             State((:θ, :Ω), initial_condition(isolated_blob, domain), domain)

        prob = SpectralODEProblem(Linear!, NonLinear!, ic, domain, tspan;
            p=parameters, dt=dt, operators=:all)

        filename = joinpath(@__DIR__, "output", "checkpoint_$(wrapper)_$(typeof(scheme)).h5")
        output = Output(prob; filename=filename, simulation_name="checkpoint_test",
            store_locally=false)

        cache = get_cache(prob, scheme)
        perform_step!(cache, prob, first(tspan))

        save_checkpoint!(output, cache, 1, first(tspan) + dt)
        restored = restore_checkpoint(output.simulation, prob, scheme)

        @test typeof(restored.u) == typeof(cache.u)
        @test get_domain(restored.u) === domain
        @test get_data(restored.u) ≈ get_data(cache.u)

        close(output.simulation.file)
        rm(filename; force=true)
    end
end
