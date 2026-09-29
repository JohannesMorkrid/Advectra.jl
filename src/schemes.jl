# TODO inherit type from SciML 
abstract type AbstractODEAlgorithm end

# Equivalent to OrdinaryDiffEqConstantCache in SciML
abstract type AbstractTableau end

# Equivalent to DECache in SciMLBase 
abstract type AbstractCache end

function get_cache(prob::SpectralODEProblem, alg::AbstractODEAlgorithm)
    get_cache(prob, alg, isinplace(prob))
end

# ----------------------------------------- MSS1 -------------------------------------------

struct MSS1 <: AbstractODEAlgorithm end

mutable struct MSS1ConstantCache{U, C} <: AbstractCache
    #Coefficents are all 1
    u::U
    c::C
end

# Constant cache
function get_cache(prob::SpectralODEProblem, alg::MSS1, ::Val{false})
    dt = prob.dt
    u = copy(prob.u0_hat)
    # Calculate linear differential operator coefficent once to cache it
    D = prob.L(one.(similar(prob.u0_hat)), prob.p, 0)
    c = @. (1 - D * dt)^-1
    MSS1ConstantCache(u, c)
end

@muladd function perform_step!(cache::MSS1ConstantCache, prob::SpectralODEProblem,
                               t::Number)
    @unpack u, c = cache
    @unpack N, p, dt = prob

    # TODO perhaps relax assumption on D being a constant coefficient

    # Perform step (.* needed for elementwise multiplication)
    cache.u = c .* (u + dt * N(u, p, t))
end

mutable struct MSS1Cache{U, C, K} <: AbstractCache
    #Coefficents are all 1
    u::U
    u_prev::U
    c::C
    k::K
end

# In-place cache
function get_cache(prob::SpectralODEProblem, alg::MSS1, ::Val{true})
    dt = prob.dt
    u = copy(prob.u0_hat)
    u_prev = copy(u)
    D = similar(u)
    # Calculate linear differential operator coefficent once to cache it
    prob.L(D, one.(similar(prob.u0_hat)), prob.p, 0)
    c = @. (1 - D * dt)^-1
    k = zero(D)
    MSS1Cache(u, u_prev, c, k)
end

@muladd function perform_step!(cache::MSS1Cache, prob::SpectralODEProblem, t::Number)
    @unpack u, u_prev, c, k = cache
    @unpack N, p, dt = prob

    # Compute difference
    N(k, u_prev, p, t)
    # Perform step in-place
    @. u = c * (u_prev + dt * k) # TODO fix lock conflict here!
    u_prev .= u
end

# ----------------------------------------- MSS2 -------------------------------------------

struct MSS2 <: AbstractODEAlgorithm end

# Coefficents/tableu
struct MSS2Tableau{T} <: AbstractTableau
    g0::T # = 3 / 2
    a0::T # = -1 / 2
    a1::T # = 2
    b0::T # = -1
    b1::T # = 2
end

function MSS2Tableau{T}() where {T}
    g0 = 3 // 2
    a0 = -1 // 2
    a1 = 2
    b0 = -1
    b1 = 2
    MSS2Tableau{T}(g0, a0, a1, b0, b1)
end

MSS2Tableau() = MSS2Tableau{Float64}()

mutable struct MSS2ConstantCache{U, C, K} <: AbstractCache
    u::U
    c::C
    u0::U
    k0::K
    tab::MSS2Tableau
    step::Int
end

function get_cache(prob::SpectralODEProblem, alg::MSS2, ::Val{false})
    tab = MSS2Tableau{get_precision(prob)}()
    dt = prob.dt
    u = copy(prob.u0_hat)
    D = prob.L(one.(similar(prob.u0_hat)), prob.p, 0)
    c = @. (tab.g0 - D * dt)^-1
    u0 = copy(prob.u0_hat)
    k0 = zero(u)
    step = 1
    MSS2ConstantCache(u, c, u0, k0, tab, step)
end

@muladd function perform_step!(cache::MSS2ConstantCache, prob::SpectralODEProblem,
                               t::Number)
    @unpack u, c, u0, k0, tab = cache
    @unpack a0, a1, b0, b1 = tab
    @unpack N, p, dt = prob

    if cache.step == 1
        cache.step += 1
        # Perform step using MSS1
        cache1 = get_cache(prob, MSS1())
        perform_step!(cache1, prob, t)
        cache.k0 = N(u0, p, t)
        cache.u = cache1.u
    else
        k1 = N(u, p, t)
        # Step
        cache.u = c .* (a0 * u0 + a1 * u + dt * (b0 * k0 + b1 * k1))
        # Shifting values downwards    
        cache.k0 = k1
        cache.u0 = u
    end
end

mutable struct MSS2Cache{U, C, K} <: AbstractCache
    u::U
    c::C
    u0::U
    u1::U
    k0::K
    k1::K
    tab::MSS2Tableau
    step::Int
end

function get_cache(prob::SpectralODEProblem, alg::MSS2, ::Val{true})
    tab = MSS2Tableau{get_precision(prob)}()
    dt = prob.dt
    u = copy(prob.u0_hat)
    D = similar(u)
    prob.L(D, one.(similar(prob.u0_hat)), prob.p, 0)
    c = @. (3 / 2 - D * dt)^-1
    u0 = copy(prob.u0_hat)
    u1 = zero(u)
    k0 = zero(u)
    k1 = zero(u)
    step = 1
    MSS2Cache(u, c, u0, u1, k0, k1, tab, step)
end

@muladd function perform_step!(cache::MSS2Cache, prob::SpectralODEProblem, t::Number)
    @unpack u, c, u0, u1, k0, k1, tab = cache
    @unpack a0, a1, b0, b1 = tab
    @unpack N, p, dt = prob

    if cache.step == 1
        cache.step += 1
        # Perform step using MSS1
        cache1 = get_cache(prob, MSS1())
        perform_step!(cache1, prob, t)
        N(k0, u0, p, t)
        u1 .= cache1.u
        u .= u1
    else
        N(k1, u1, p, t)
        # Step
        @. u = c * (a0 * u0 + a1 * u1 + dt * (b0 * k0 + b1 * k1))
        # Shifting values downwards    
        k0 .= k1
        u0 .= u1
        u1 .= u
    end
end

# ----------------------------------------- MSS3 -------------------------------------------

struct MSS3 <: AbstractODEAlgorithm end

# Coefficents/tableu
struct MSS3Tableau{T} <: AbstractTableau
    g0::T # = 11 / 6
    a0::T # = 1 / 3
    a1::T # = -3/2
    a2::T # = 3
    b0::T # = 1
    b1::T # = -3
    b2::T # = 3
end

function MSS3Tableau{T}() where {T}
    g0 = 11 // 6
    a0 = 1 // 3
    a1 = -3 // 2
    a2 = 3
    b0 = 1
    b1 = -3
    b2 = 3
    MSS3Tableau{T}(g0, a0, a1, a2, b0, b1, b2)
end

MSS3Tableau() = MSS3Tableau{Float64}()

mutable struct MSS3ConstantCache{U, C, K} <: AbstractCache
    u::U
    c::C
    u0::U
    u1::U
    k0::K
    k1::K
    tab::MSS3Tableau
    step::Int
end

function get_cache(prob::SpectralODEProblem, alg::MSS3, ::Val{false})
    tab = MSS3Tableau{get_precision(prob)}()
    dt = prob.dt
    u = copy(prob.u0_hat)
    D = prob.L(one.(similar(prob.u0_hat)), prob.p, 0)
    c = @. (tab.g0 - D * dt)^-1
    u0 = copy(prob.u0_hat)
    u1 = zero(u)
    k0 = zero(u)
    k1 = zero(u)
    step = 1
    MSS3ConstantCache(u, c, u0, u1, k0, k1, tab, step)
end

@muladd function perform_step!(cache::MSS3ConstantCache, prob::SpectralODEProblem,
                               t::Number)
    @unpack u, c, u0, u1, k0, k1, tab = cache
    @unpack a0, a1, a2, b0, b1, b2 = tab
    @unpack N, p, dt = prob

    if cache.step == 1
        cache.step += 1
        # Perform step using MSS1
        cache1 = get_cache(prob, MSS1())
        perform_step!(cache1, prob, t)
        cache.k0 = N(u0, p, t)
        cache.u1 = cache1.u
        cache.u = cache.u1 # For output handling
    elseif cache.step == 2
        cache.step += 1
        # Set up cache for stepping with MSS2
        cache2 = get_cache(prob, MSS2())
        cache2.step = 2
        cache2.k0 = k0
        cache2.u = u1
        # Perform step using MSS2
        perform_step!(cache2, prob, t)
        cache.k1 = cache2.k0
        cache.u = cache2.u
    else
        k2 = N(u, p, t)
        # Step
        cache.u = c .* (a0 * u0 + a1 * u1 + a2 * u + dt * (b0 * k0 + b1 * k1 + b2 * k2))
        # Shifting values downwards    
        cache.k0 = k1
        cache.k1 = k2
        cache.u0 = u1
        cache.u1 = u
    end
end

mutable struct MSS3Cache{U, C, K} <: AbstractCache
    u::U
    c::C
    u0::U
    u1::U
    u2::U
    k0::K
    k1::K
    k2::K
    tab::MSS3Tableau
    step::Int
end

function get_cache(prob::SpectralODEProblem, alg::MSS3, ::Val{true})
    tab = MSS3Tableau{get_precision(prob)}()
    dt = prob.dt
    u = copy(prob.u0_hat)
    D = similar(u)
    prob.L(D, one.(similar(prob.u0_hat)), prob.p, 0)
    c = @. (tab.g0 - D * dt)^-1
    u0 = copy(prob.u0_hat)
    u1 = zero(u)
    u2 = zero(u)
    k0 = zero(u)
    k1 = zero(u)
    k2 = zero(u)
    step = 1
    MSS3Cache(u, c, u0, u1, u2, k0, k1, k2, tab, step)
end

@muladd function perform_step!(cache::MSS3Cache, prob::SpectralODEProblem, t::Number)
    @unpack u, c, u0, u1, u2, k0, k1, k2, tab = cache
    @unpack a0, a1, a2, b0, b1, b2 = tab
    @unpack N, p, dt = prob

    if cache.step == 1
        cache.step += 1
        # Perform step using MSS1
        cache1 = get_cache(prob, MSS1())
        perform_step!(cache1, prob, t)
        N(k0, u0, p, t)
        u1 .= cache1.u
        u .= u1 # For output handling
    elseif cache.step == 2
        cache.step += 1
        # Set up cache for stepping with MSS2
        cache2 = get_cache(prob, MSS2())
        cache2.step = 2
        cache2.k0 .= k0
        cache2.u1 .= u1
        # Perform step using MSS2
        perform_step!(cache2, prob, t)
        k1 .= cache2.k0
        u2 .= cache2.u
        u .= u2 # For output handling
    else
        N(k2, u2, p, t)
        # Step
        @. u = c * (a0 * u0 + a1 * u1 + a2 * u2 + dt * (b0 * k0 + b1 * k1 + b2 * k2))
        # Shifting values downwards    
        k0 .= k1
        k1 .= k2
        u0 .= u1
        u1 .= u2
        u2 .= u
    end
end

# TODO investigate if possible to remove u1 for MSS2 and similarly u2 for MSS3
# TODO perhaps fix return type?