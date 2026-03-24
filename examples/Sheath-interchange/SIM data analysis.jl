## Run all (alt+enter)
using Advectra

domain = Domain(256, 256, 48, 48; dealiased=true)
using Statistics
# Open data file
fid = h5open("output/sheath-interchange long time series.h5", "r")

## Open solution
S = 1
for S in eachindex(keys(fid))
    try
        simulation = fid[keys(fid)[S]]
        try
            data = read(simulation["Density probe/data"])
            t = read(simulation["Density probe/t"])
            display(plot(t, data'; xlabel=L"t", ylabel=L"n(0,0)"))
        catch
        end
    catch
    end
end

close(fid)

fid = h5open("output/sheath-interchange april forth.h5", "r")
simulation = fid[keys(fid)[1]]
read_attribute(simulation, "g")

n = read(simulation["Density probe/data"])
n = stack(n)

using Statistics
using StatsPlots
using Distributions

n_n = (n .- mean(n)) / std(n)

density(n_n[1:500001]; minorticks=true,
        xlabel=L"(\tilde{n}-\langle\tilde{n}\rangle)\tilde{n}_{rms}",
        ylabel=L"P(\tilde{n})", label="", frame=:box, figsize=(3.37, 2.08277),
        fontfamily="Computer Modern",
        linewidth=0.75, grid=false)
plot!(Normal(0, 1); label=L"N(0,1)", guidefontsize=15, tickfontsize=15, legendfontsize=13,
      markersize=2.25, linewidth=0.75, minorticks=4)
savefig("fluctuation pdf.pdf")

density(n)

keys(fid)
S = 1
simulation = fid[keys(fid)[S]]
data = read(simulation["Density probe/data"])
t = read(simulation["Density probe/t"])
display(plot(t, data'; xlabel=L"t", ylabel=L"n(0,0)"))
#delete_object(simulation)

data = simulation["fields"][:, :, :, :]
t = simulation["t"][:]
default(; legend=false)
anim = @animate for i in axes(data, 4)
    heatmap(data[:, :, 1, i]; aspect_ratio=:equal, xaxis=L"x", yaxis=L"y",
            title=L"n(t=" * "$(round(t[i], digits=0)))")
end
gif(anim, "delete me.gif"; fps=20)

simulation

using Plots
default(; legend=true, tickfont=font(12))
x = -6:0.1:6
plots = map(f -> plot(x, f.(x); title=string(f), tickfont=font(12)), [sin, cos, sinh, cosh])
plotgrid = plot(plots...; layout=grid(2, 2), link=:both)

#
using JLD
save("density probe.jld", "probe data", n)

# Extract data to do local python analysis
using Advectra
fid = h5open("output/sheath-interchange g=1e-3.h5", "r")

data = fid[keys(fid)[1]]["All probe/data"][:, :, 1:1:end]
t = fid[keys(fid)[1]]["All probe/t"][1:10:end]

using JLD
jldopen("output/all probes g=5e-3 dense.jld", "w") do file
    g = create_group(file, "data")
    g["n"] = data[1, 1, 1:argmax(t)]
    g["Omega"] = data[1, 2, 1:argmax(t)]
    g["phi"] = data[1, 3, 1:argmax(t)]
    g["vx"] = data[1, 4, 1:argmax(t)]
    g["Gamma"] = data[1, 5, 1:argmax(t)]
    g["t"] = t[1:argmax(t)]
end

fields = fid[keys(fid)[1]]["fields"][:, :, :, :]
heatmap(fields[:, :, 1, 25]; aspect_ratio=:equal)

@views n = data[1, 1, :]

plot(n[500:10:end]; marker=".")

data = fid[keys(fid)[1]]["Radial flux/data"][:, :, :]
data = fid[keys(fid)[1]]["Enstropy energy integral/data"][:]
data = fid[keys(fid)[1]]["Potential energy integral/data"][:]
data = fid[keys(fid)[1]]["All probe/data"][:, :, 1:10:end]

"Radial flux"
"Kinetic energy integral"

#data = fid[keys(fid)[1]]["Radial flux/data"][1:10:5000001]
data = fid[keys(fid)[1]]["Enstropy energy integral/data"][:]

plot(data[(end-100):1:end]; marker=".")
fid[keys(fid)[1]]

data = fid[keys(fid)[1]]["fields"][:, :, :, :]
heatmap(data[:, :, 1, end])#,aspect_ratio=:equal)

using Advectra
fid = h5open("output/gyro-bohm=1e-2 CUDA.h5", "r")
sim = fid[keys(fid)[1]]
data = sim["All probe/data"][:, :, 1:N]
t = sim["All probe/t"][1:N]
dt = 0.01
t = LinRange(0.0, dt * N, N)

sim["All probe/t"][N+12_987_702]
N
N = sim["t"][2] * argmax(sim["t"][:]) / sim["All probe/t"][2]

N = round(Int64, N)

using JLD
jldopen("output/all probes sigma=1e-2 10 probes CUDA.jld", "w") do file
    g = create_group(file, "data")
    g["n"] = data[:, 1, 1:argmax(t)]
    g["Omega"] = data[:, 2, 1:argmax(t)]
    g["phi"] = data[:, 3, 1:argmax(t)]
    g["vx"] = data[:, 4, 1:argmax(t)]
    g["Gamma"] = data[:, 5, 1:argmax(t)]
    g["t"] = t[1:argmax(t)]
end

## Thesis images
fid = h5open("output/gyro-bohm=1e-1 CUDA.h5", "r")
domain = Domain(256, 256, 48, 48; dealiased=true, use_cuda=false)

points = [(x, 0) for x in range(-24, 19.2, 10)]

for N in 30
    sim = fid[keys(fid)[1]]
    sigma = read_attribute(sim, "sigma")
    n = sim["fields"][:, :, 1, N]
    t = sim["t"][N]
    heatmap(domain, n; title="Probe positions", xlabel=L"x\ [\rho_s]",
            ylabel=L"y\ [\rho_s]", size=[600, 550], margin=0Plots.px,
            top_margin=-100Plots.px, bottom_margin=-40Plots.px, titlefontsize=12,
            labelfontsize=10)#, color="black")
    #display(contour!(domain.x, domain.y, n, color=:black))
    scatter!(points; label="", markersize=4.5)
    display(plot!())
    savefig("probe positions.pdf")
end

sim = fid[keys(fid)[1]]

K = sim["Kinetic energy integral/data"][1:320_959]
P = sim["Potential energy integral/data"][1:320_959]

probe_data = sim["All probe/data"][:, :, 1:320_959]
n = probe_data[:, 1, :]
Ω = probe_data[:, 2, :]
ϕ = probe_data[:, 3, :]
vx = probe_data[:, 4, :]
Γ = probe_data[:, 5, :]

plot(K; aspect_ratio=:auto)
plot!(P)

using Statistics
K_n = (K .- mean(K)) ./ (std(K))
P_n = (P .- mean(P)) ./ std(P)
n_n = (n .- mean(n; dims=1)) ./ std(n; dims=1)
ϕ_n = (ϕ .- mean(ϕ; dims=1)) ./ std(ϕ; dims=1)
Ω_n = (Ω .- mean(Ω; dims=1)) ./ std(Ω; dims=1)
vx_n = (vx .- mean(vx; dims=1)) ./ std(vx; dims=1)
Γ_n = (Γ .- mean(Γ; dims=1)) ./ std(Γ; dims=1)

plot(K_n; aspect_ratio=:auto)
plot!(P_n; aspect_ratio=:auto)

plot(K_n[1:10000]; aspect_ratio=:auto)
plot(P_n[1:10000]; aspect_ratio=:auto)

plot(n_n[:, 1:7000]'; aspect_ratio=:auto)
plot(ϕ_n[:, 1:10000]'; aspect_ratio=:auto)
plot(Ω_n[:, 1:10000]'; aspect_ratio=:auto)

histogram(n_n'; nbis=64, aspect_ratio=:auto, yaxis=:log10)
histogram(ϕ_n'; nbis=64, aspect_ratio=:auto, yaxis=:log10)
histogram(Ω_n'; nbis=64, aspect_ratio=:auto, yaxis=:log10)
histogram(vx_n'; nbis=64, aspect_ratio=:auto, yaxis=:log10)
histogram(Γ_n'; nbis=64, aspect_ratio=:auto)

sim["All probe/t"][7000]

spectra = sim["Radial potential energy spectra/data"][:, 1, 1:100_000]
sim["Radial kinetic energy spectra/t"][10000]

mean_spectra = mean(spectra; dims=2)

plot(mean_spectra; xaxis=:log10, aspect_ratio=:auto, label="")

heatmap(sim["fields"][:, :, 1, 100])
modes = fft(sim["fields"][:, :, 1, 100])

heatmap(sign.(imag.(modes)))

sign(-2)

plot(fftshift(domain.kx)[130:end], fftshift(mean_spectra)[130:end]; aspect_ratio=:auto,
     xaxis=:log10)

fftshift(domain.kx)[129]
argmax(fftshift(mean_spectra))

sim

heatmap(fftshift(abs.(modes)))

plot(abs.(modes)[:, 1]; aspect_ratio=:auto, yaxis=:log10)

## Determine structure size
using Statistics
n = sim["fields"][:, :, 1, :]
n_x_2 = sum(n; dims=1)[1, :, :] / domain.Ny
specter1 = mean([abs.(fft(n_x_1[:, i])) for i in 1:1001])
specter2 = mean([abs.(fft(n_x_2[:, i])) for i in 1:1001])

plot(fftshift(specter1)[129:end] / sum(fftshift(specter1)[129:end]); aspect_ratio=:auto)
plot!(fftshift(specter2)[129:end] / sum(fftshift(specter2)[129:end]); aspect_ratio=:auto)

heatmap()

plot(abs.(fft(mean(n; dims=3)[:, :, 1])); aspect_ratio=:auto, label="")
fid = h5open("output/sheath-interchange g=1e-2.h5", "r")
sim = fid["D_n=0.01, D_Ω=0.01, N=1.0, g=0.01, kappa=0.31622776601683794, sigma_n=0.001, sigma_Ω=0.001"]

## ------------------------------- Gyro Bohm -----------------------------------------------
fid = h5open("output/gyro-bohm=5e-2.h5", "r")
sim = fid["D=0.01, g=0.1, sigma=0.05"]
t = sim["All probe/t"][1:3_079_300]
data = sim["All probe/data"][:, :, 1:3_079_300]

sim["All probe/t"][3_079_300]

argmax(t)
heatmap(sim["fields"][:, :, 1, argmax(sim["t"][:])]; levels=17)#, cmap=:black)

using JLD
jldopen("output/all probes gyro-bohm=5e-2 10 probes.jld", "w") do file
    g = create_group(file, "data")
    g["n"] = data[:, 1, 1:argmax(t)]
    g["Omega"] = data[:, 2, 1:argmax(t)]
    g["phi"] = data[:, 3, 1:argmax(t)]
    g["vx"] = data[:, 4, 1:argmax(t)]
    g["Gamma"] = data[:, 5, 1:argmax(t)]
    g["t"] = t[1:argmax(t)]
end
