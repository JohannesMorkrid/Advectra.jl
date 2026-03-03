module AdvectraComponentArraysExt

using Advectra, ComponentArrays, FFTW

import Advectra: _allocate_coefficients
function _allocate_coefficients(u0::ComponentArray, domain::Domain)
    ComponentArray(;
                   (key => _allocate_coefficients(getproperty(u0, key), domain)
                    for key in keys(u0))...)
end

import Advectra: _spectral_transform!
function _spectral_transform!(du, u::ComponentArray, p::P) where {P<:FFTW.Plan}
    for k in keys(u)
        _spectral_transform!(getproperty(du, k), getproperty(u, k), p)
    end
end

import Advectra: assert_no_nan
assert_no_nan(u::ComponentArray, t) = assert_no_nan(parent(u), t)

# TODO perhaps custom write_state? So that output is easier to read

end