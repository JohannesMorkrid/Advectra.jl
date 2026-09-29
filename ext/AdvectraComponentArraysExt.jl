module AdvectraComponentArraysExt

using Advectra, ComponentArrays, FFTW

import Advectra: prepare_initial_condition
# Convert the underlying data only, converting the whole ComponentArray drops its axes
function prepare_initial_condition(u0::ComponentArray, domain::Domain)
    ComponentArray(getdata(u0) |> memory_type(domain, Physical()), getaxes(u0))
end

import Advectra: _allocate_coefficients
function _allocate_coefficients(u0::ComponentArray, domain::Domain)
    ComponentArray(;
                   (key => _allocate_coefficients(getproperty(u0, key), domain)
                    for key in keys(u0))...)
end

import Advectra: _spectral_transform!
function _spectral_transform!(du, p::P, u::ComponentArray) where {P<:FFTW.Plan}
    for k in keys(u)
        _spectral_transform!(getproperty(du, k), p, getproperty(u, k))
    end
end

import Advectra: assert_no_nan
assert_no_nan(u::ComponentArray, t) = assert_no_nan(parent(u), t)

# TODO perhaps custom write_state? So that output is easier to read

end