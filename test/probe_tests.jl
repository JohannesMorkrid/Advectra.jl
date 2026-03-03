# ------------------------------------------------------------------------------------------
#                                        Probe Test                                         
# ------------------------------------------------------------------------------------------

using Advectra
using CUDA
import Advectra: build_diagnostic

domain = Domain(256, 256; MemoryType=CuArray)
ic = initial_condition(isolated_blob, domain) |> Advectra.memory_type(domain)

probe = build_diagnostic(Val(:probe_all); domain=domain,
                         positions=[(0, 0), (0.1, 0), (0.4, 0)])

ic_hat = cat(get_fwd(domain) * ic[:, :, 1], get_fwd(domain) * ic[:, :, 2]; dims=3)

probe(ic_hat, (; domain=domain, operators=()), 0.0)

"""
    Test the following:
    * Does construction of probe throw an error when wrong Tuple length
    * Does construction of probe throw an error when point is outside of domain bounds.
    * Does construction of probe promote to the right type
    
    * Does interpolation=nothing lead to indicies instead of positions
    * Does it probe correctly based on:
        - GPUArray
        - Array
    * What about whether or not the Array has shape (256,256) or (256,256,1)
    * Perhaps also one test using interpolation

    * Check all 5 probe types
"""