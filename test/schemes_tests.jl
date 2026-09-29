# ------------------------------------------------------------------------------------------
#                                   Time Integration Tests
# ------------------------------------------------------------------------------------------

using Test
using Advectra
import Advectra: get_cache, perform_step!

# du/dt = ν∇²u + au, with the diffusion treated implicitly and the growth explicitly. For a
# single mode with k² = 5 the exact solution is u(t) = u(0)exp((a - 5ν)t)
const domain = Domain(8; L=2π)
const u0 = @. cos(domain.x') * sin(2 * domain.y)
const parameters = (ν=0.2, a=0.5)
exact(t) = u0 .* exp((parameters.a - 5parameters.ν) * t)

Linear!(du, u, operators, p, t) = du .= p.ν .* operators.laplacian(u)
NonLinear!(du, u, operators, p, t) = du .= p.a .* u
Linear(u, operators, p, t) = p.ν .* operators.laplacian(u)
NonLinear(u, operators, p, t) = p.a .* u

function problem(dt, inplace)
    L, N = inplace ? (Linear!, NonLinear!) : (Linear, NonLinear)
    SpectralODEProblem(L, N, u0, domain, [0.0, 1.0]; p=parameters, dt=dt)
end

# Maximum error at t = 1
function solve_error(scheme, dt, inplace)
    prob = problem(dt, inplace)
    cache = get_cache(prob, scheme)
    for step in 0:(round(Int, 1 / dt) - 1)
        perform_step!(cache, prob, step * dt)
    end
    maximum(abs, bwd_plan(domain) * cache.u .- exact(1.0))
end

order(scheme, inplace) = log2(solve_error(scheme, 0.1, inplace) /
                              solve_error(scheme, 0.05, inplace))

@testset "Convergence order - inplace=$inplace" for inplace in (true, false)
    @test order(MSS1(), inplace) ≈ 1 atol = 0.1
    @test order(MSS2(), inplace) ≈ 2 atol = 0.1
    # MSS3 starts with an MSS1 and an MSS2 step, which limits it to second order (BUGS.md)
    @test order(MSS3(), inplace) ≈ 2 atol = 0.1
end

@testset "MSS3 is third order with an exact start" begin
    function error_exact_start(dt)
        prob = problem(dt, true)
        cache = get_cache(prob, MSS3())
        # Replace the MSS1/MSS2 start-up steps with the exact history
        F(t) = fwd_plan(domain) * exact(t)
        cache.u0 .= F(0)
        cache.u1 .= F(dt)
        cache.u2 .= F(2dt)
        cache.u .= F(2dt)
        prob.N(cache.k0, cache.u0, prob.p, 0.0)
        prob.N(cache.k1, cache.u1, prob.p, dt)
        cache.step = 3
        for step in 2:(round(Int, 1 / dt) - 1)
            perform_step!(cache, prob, step * dt)
        end
        maximum(abs, bwd_plan(domain) * cache.u .- exact(1.0))
    end
    @test log2(error_exact_start(0.1) / error_exact_start(0.05)) ≈ 3 atol = 0.2
end

@testset "In-place and out-of-place agree - $(nameof(typeof(scheme)))" for
    scheme in (MSS1(), MSS2(), MSS3())

    @test solve_error(scheme, 0.1, true) ≈ solve_error(scheme, 0.1, false)
end

@testset "Float32 precision" begin
    # The same run in single and double precision
    solutions = map((Float32, Float64)) do T
        d = Domain(8; L=2π, precision=T)
        prob = SpectralODEProblem(Linear!, NonLinear!, @.(cos(d.x') * sin(2 * d.y)), d,
                                  [0.0, 1.0]; p=parameters, dt=0.05)
        cache = get_cache(prob, MSS3())
        for step in 0:19
            perform_step!(cache, prob, step * 0.05)
        end
        cache.u
    end
    @test eltype(first(solutions)) == ComplexF32
    @test first(solutions) ≈ last(solutions) rtol = 1e-5
end
