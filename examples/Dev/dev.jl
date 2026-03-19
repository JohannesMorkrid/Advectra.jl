using Advectra
using CUDA

# ------------------------------------------------------------------------------------------
#                                      Specify Domain                                       
# ------------------------------------------------------------------------------------------

domain = Domain(256, 256; Lx=50, Ly=50, MemoryType=CuArray, precision=Float64)

# ------------------------------------------------------------------------------------------
#                                       Specify RHSs                                        
# ------------------------------------------------------------------------------------------

Linear(du, u, model, t) = "..."

function NonLinear(du, u, model, t)
    @context model
    u = 0
end

# ------------------------------------------------------------------------------------------
#                                     Clarify Variables                                     
# ------------------------------------------------------------------------------------------

variables = (vx=:spectral, Vx=:physical, n=:density, Ω=:vorticity, ϕ=:potential)

ic = (n=:blob, Ω=:zero)

initial_condition!(variables.n, :blob; "...")
initial_condition!(variables.Ω, :blob; "...")

n = initial_conditon()
Ω = initial_condition()
ϕ = :potential

variables = (n=:blob → :density,
             Ω=:zeros → :vorticity,
             ϕ=:spectral → :potential)

Variable(; symbol=:n, data=:pointer, context=:density)

variables.n => Variable(:n, data, :density)

Field(:n, data → u[:, :, 1], :density, domain)
Field(:n, data, :density, domain)

# Programmatic way:
struct VariablesHardCoded{physical,spectral}
    n::spectral
    N::physical
    ϕ::spectral

    function Variables(domain)
        n = initial_condition(blob, domain; A=1, B=1)
        N = zero(domain)
        ϕ = spectral_zero(domain)
    end
end

[
    Variable(:n, initial_condition(blob, domain); context=:density),
    Variable(:Ω, initial_condition(blob, domain); context=:vorticity),
    Variable(:ϕ, spectral_zero)
]

(n, Ω, density, vorticity)

ic = [:n, :Omega, :ϕ]

# ------------------------------------------------------------------------------------------
#                              Specify Parameters And Timespan                              
# ------------------------------------------------------------------------------------------

parameters = (ν=1e-2, κ=1e-2, A=A)

# Time interval
tspan = [0.0, 20.0]

# ------------------------------------------------------------------------------------------
#                                    Specify Diagnostics                                    
# ------------------------------------------------------------------------------------------

Diagnostic(probe, :all; stride=100, positions=[(5, 0), (8.5, 0), (11.25, 0), (14.375, 0)])
Diagnostic(probe, :density; stride=100,
           positions=[(5, 0), (8.5, 0), (11.25, 0), (14.375, 0)])
Diagnostic(probe, :n; stride=100, positions=[(5, 0), (8.5, 0), (11.25, 0), (14.375, 0)])

probe(:density)
sample()
spectra()
display()

method::M

# Can be determined from method
name::N
metadata::L
assumes_spectral_state::Bool
stores_data::Bool
args::A
kwargs::K

→(A, B) = 1
# Array of diagnostics want
diagnostics = [
    :sample → (; stride=1000, fields=[:density, :vorticity, :potential]),
    :plot → (; stride=1000, fields=[:density, :vorticity, :potential]),
    :probe →
    (; positions=[(5, 0), (8.5, 0), (11.25, 0), (14.375, 0)], stride=10, fields=:all),
    :progress → (; stride=-1),
    :cfl → (; stride=250, silent=true, storage_limit="2KB")
]

# ------------------------------------------------------------------------------------------
#                               Construct SpectralODEProblem                                
# ------------------------------------------------------------------------------------------

prob = SpectralODEProblem(Linear, NonLinear, ic, domain, tspan; p=parameters, dt=2.5e-3,
                          boussinesq=false, diagnostics=diagnostics, density=:linear,
                          operators=:all)

# ------------------------------------------------------------------------------------------
#                                     Construct Output                                      
# ------------------------------------------------------------------------------------------

output = Output(prob; filename="", resume=true)

# ------------------------------------------------------------------------------------------
#                                      Run Simulation                                       
# ------------------------------------------------------------------------------------------

sol = solve(::SpectralODEProblem, MSS3(), output)