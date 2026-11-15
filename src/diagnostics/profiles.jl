# ------------------------------------------------------------------------------------------
#                                         Profiles                                          
# ------------------------------------------------------------------------------------------

#1D profile: n_0(x,t) = 1/L_y∫_0^L_y n(x,y,t)dy
function radial_density_profile(state, prob, time; quadrature=nothing)
    sum(selectdim(state, ndims(prob.domain), 1); dims=1)' ./ prob.domain.Ly
end

function poloidal_density_profile(state, prob, time; quadrature=nothing)
    sum(selectdim(state, ndims(prob.domain), 1); dims=2) ./ prob.domain.Lx
end

function radial_vorticity_profile(state, prob, time; quadrature=nothing)
    sum(selectdim(state, ndims(prob.domain), 2); dims=1)' ./ prob.domain.Ly
end

function poloidal_vorticity_profile(state, prob, time; quadrature=nothing)
    sum(selectdim(state, ndims(prob.domain), 2); dims=2) ./ prob.domain.Lx
end

#Γ_0(x,t) = 1/L_y∫_0^L_y nv_x dy
function radial_flux_profile(state, prob, time; quadrature=nothing)
    v_x, v_y = vExB(state, prob.domain)
    sum(selectdim(state, ndims(prob.domain), 1) .* v_x; dims=1)' / prob.domain.Ly
end
