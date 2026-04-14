# ------------------------------------------------------------------------------------------
#                                       QuadraticTerm                                       
# ------------------------------------------------------------------------------------------

"""
"""
struct QuadraticTerm{TP<:AbstractTransformPlans,P<:AbstractArray,SP<:AbstractArray,
                     C<:AbstractFloat} <: NonLinearOperator
    transforms::TP
    U::P
    V::P
    up::SP
    vp::SP
    padded::Bool
    dealiasing_coefficient::C

    function QuadraticTerm(domain::AbstractDomain)
        precision = get_precision(domain)
        Ns = size(domain)

        # Allocate the physical array, using zero-padding if dealiasing is enabled
        utmp = zeros(precision, domain.dealiased ? pad_size(Ns) : Ns) |> memory_type(domain)

        transforms = construct_transform_plans(utmp, Domain, Val(domain.real_transform))

        # Allocate data for pseudo spectral schemes
        up = fwd(transforms) * zero(utmp)
        vp = similar(up)
        U = bwd(transforms) * up
        V = similar(U)

        # Calculate correct conversion coefficent
        dealiasing_coefficient = precision(length(up) / spectral_length(domain))

        new{typeof(transforms),typeof(U),typeof(up),
            typeof(dealiasing_coefficient)}(transforms, U, V, up, vp,
                                            domain.dealiased, dealiasing_coefficient)
    end
end

"""
    pad_size(Ns::NTuple{N,Int})
  
  Computes the zero-pad size for the size Tuple ``Ns`` following the 3/2 rule.
"""
function pad_size(Ns::NTuple{N,Int}) where {N}
    ntuple(i -> Ns[i] > 1 ? div(3 * Ns[i], 2, RoundUp) : 1, N)
end
# -------------------------------- Constructor related -------------------------------------

build_operator(::Val{:quadratic_term}, domain::Domain; kwargs...) = QuadraticTerm(domain)

#------------------------------ Quadratic terms interface ----------------------------------

# In-place operator
@inline function (quadratic_term::QuadraticTerm)(out::T, u::T,
                                                 v::T) where {T<:AbstractGPUArray}
    @assert size(u)==size(v) "u and v must have the same size"
    @unpack transforms, U, V, up, vp, padded, dealiasing_coefficient = quadratic_term

    mul!(U, bwd(transforms), padded ? pad!(up, u, typeof(transforms)) : u)
    mul!(V, bwd(transforms), padded ? pad!(vp, v, typeof(transforms)) : v)
    @. U *= V

    mul!(padded ? up : out, fwd(transforms), U)
    padded ? out .= dealiasing_coefficient .* unpad!(out, up, typeof(transforms)) : out
end

@inline function (quadratic_term::QuadraticTerm)(out::AbstractArray, u::AbstractArray,
                                                 v::AbstractArray)
    @assert size(u)==size(v) "u and v must have the same size"
    @unpack transforms, U, V, up, vp, padded, dealiasing_coefficient = quadratic_term

    # Spawn threads to perform mul! in parallel
    task_U = Threads.@spawn mul!(U, bwd(transforms),
                                 padded ? pad!(up, u, typeof(transforms)) : u)
    task_V = Threads.@spawn mul!(V, bwd(transforms),
                                 padded ? pad!(vp, v, typeof(transforms)) : v)
    # Wait for both tasks to finish
    wait(task_V)
    wait(task_U)

    @threads for i in eachindex(U)
        U[i] *= V[i]
    end
    mul!(padded ? up : out, fwd(transforms), U)
    padded ? out .= dealiasing_coefficient .* unpad!(out, up, typeof(transforms)) : out
end

# Out-of-place operator
(op::QuadraticTerm)(u::AbstractArray, v::AbstractArray) = op(similar(u), u, v)

# ------------------------------------------ Helpers ---------------------------------------

# Specialized for 2D arrays
@inline function pad!(up::DU, u::U,
                      ::Type{<:FFTPlans}) where {T,DU<:AbstractArray{T},U<:AbstractArray{T}}
    up .= zero.(up)
    Ny, Nx = size(u)

    Nxl = div(Nx, 2, RoundUp)
    Nxu = div(Nx, 2, RoundDown)
    Nyl = div(Ny, 2, RoundUp)
    Nyu = div(Ny, 2, RoundDown)

    @views @inbounds up[1:Nyl, 1:Nxl] .= u[1:Nyl, 1:Nxl] # Lower left
    @views @inbounds up[1:Nyl, (end-Nxu+1):end] .= u[1:Nyl, (end-Nxu+1):end] # Lower right
    @views @inbounds up[(end-Nyu+1):end, 1:Nxl] .= u[(end-Nyu+1):end, 1:Nxl] # Upper left
    @views @inbounds up[(end-Nyu+1):end, (end-Nxu+1):end] .= u[(end-Nyu+1):end,
                                                               (end-Nxu+1):end] # Upper right
    return up
end

@inline function pad!(up::DU, u::U,
                      ::Type{<:rFFTPlans}) where {T,DU<:AbstractArray{T},
                                                  U<:AbstractArray{T}}
    up .= zero.(up)
    Ny, Nx = size(u)

    Nxl = div(Nx, 2, RoundUp)
    Nxu = div(Nx, 2, RoundDown)

    @views @inbounds up[1:Ny, 1:Nxl] .= u[1:Ny, 1:Nxl] # Lower left
    @views @inbounds up[1:Ny, (end-Nxu+1):end] .= u[1:Ny, (end-Nxu+1):end] # Lower right
    return up
end

@inline function unpad!(u::DU, up::U,
                        ::Type{<:FFTPlans}) where {T,DU<:AbstractArray{T},
                                                   U<:AbstractArray{T}}
    Ny, Nx = size(u)

    Nxl = div(Nx, 2, RoundUp)
    Nxu = div(Nx, 2, RoundDown)
    Nyl = div(Ny, 2, RoundUp)
    Nyu = div(Ny, 2, RoundDown)

    @views @inbounds u[1:Nyl, 1:Nxl] .= up[1:Nyl, 1:Nxl] # Lower left
    @views @inbounds u[1:Nyl, (end-Nxu+1):end] .= up[1:Nyl, (end-Nxu+1):end] # Lower right
    @views @inbounds u[(end-Nyu+1):end, 1:Nxl] .= up[(end-Nyu+1):end, 1:Nxl] # Upper left
    @views @inbounds u[(end-Nyu+1):end, (end-Nxu+1):end] .= up[(end-Nyu+1):end,
                                                               (end-Nxu+1):end] # Upper right
    return u
end

@inline function unpad!(u::DU, up::U,
                        ::Type{<:rFFTPlans}) where {T,DU<:AbstractArray{T},
                                                    U<:AbstractArray{T}}
    Ny, Nx = size(u)

    Nxl = div(Nx, 2, RoundUp)
    Nxu = div(Nx, 2, RoundDown)

    @views @inbounds u[1:Ny, 1:Nxl] .= up[1:Ny, 1:Nxl] # Lower left
    @views @inbounds u[1:Ny, (end-Nxu+1):end] .= up[1:Ny, (end-Nxu+1):end] # Lower right
    return u
end