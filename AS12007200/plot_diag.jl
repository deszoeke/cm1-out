# plot CM1 output

using Pkg; Pkg.activate("..")

using Revise
using NCDatasets
using Statistics
using Dates
using PyPlot

# extend PyCall masked arrays
using PyCall
using PyCall: PyObject

# PyObject method interprets Array{Union{T,Missing}} as a
# numpy masked array.
# This allows for plotting with missing values.
function PyObject(a::Array{Union{T,Missing},N}) where {T,N}
    numpy_ma = PyCall.pyimport("numpy").ma
    pycall(numpy_ma.array, Any, coalesce.(a,zero(T)), mask=ismissing.(a))
end

include("./Thermo.jl")
using .Thermo

# helper functions
pd = permutedims

# constants
const g = 9.8

filenames = filter(f -> startswith(f,"cm1out") && endswith(f,".nc"), readdir())

ds = NCDataset("data/cm1out_diag_s.nc")
zs = ds["lev"][:]
t = year.(ds["time"][:]) # hours
qv = ds["qv"][1,1,:,:]
qc = ds["qc"][1,1,:,:] # lon, lat, lev, time
prs = ds["prs"][1,1,:,:] # pressure, Pa
th = ds["th"][1,1,:,:]
the = ds["the"][1,1,:,:]
brz = ds["brz"][1,1,:,:] # bulk Richardson number as function of h
rtke = ds["rtke"][1,1,:,:]
#n2e = 9.8*diff(log.(the), dims=1)./diff(zs) # equivalent buoyancy freqency squared, using theta_e
temp = th .* (prs./1e5).^(Thermo.Rd/Thermo.Cp)

dw = NCDataset("data/cm1out_diag_w.nc")
zw = dw["lev"][:]
Ri = dw["rinum"][1,1,:,:]
wprof = dw["wprof"][1,1,:,:] # zero in test runs
nm = dw["nm"][1,1,:,:] # N^2
thinterp = dw["thinterp"][1,1,:,:]
prsinterp = dw["prsinterp"][1,1,:,:]
tempinterp = thinterp .* (prsinterp./1e5).^(Thermo.Rd/Thermo.Cp)
# n2moist = g/tempinterp * (diff(temp, dims=1)./diff(zs) - moistad(prsinterp,tempinterp)) # on w-levels
# account for latent heating effect on buoyancy:
# Nsq_moist = Nsq - g^2/(Cp*T) * gamma/(1+gamma) # avoids differencing
n2moist = nm
Ri_moist = Ri
for i in CartesianIndices((axes(nm,1)[2:end-1], axes(nm,2)))
    n2moist[i] = nm[i] + Thermo.Nsq_moist_add(tempinterp[i],prsinterp[i])
    Ri_moist[i] = (1 + Thermo.Nsq_moist_add(tempinterp[i],prsinterp[i])./nm[i]) * Ri[i]
end
# for i in 2:size(nm,1)-1
#     for j in axes(nm,2)
#         n2moist[i,j] = nm[i,j] + Thermo.Nsq_moist_add(tempinterp[i,j],prsinterp[i,j])
#         Ri_moist[i,j] = (1 + Thermo.Nsq_moist_add(tempinterp[i,j],prsinterp[i,j])./nm[i,j]) * Ri[i,j]
#     end
# end

# Ri_moist = (1 + Thermo.Nsq_moist_add.(tempinterp,prsinterp)./nm) .* Ri # factor accounts for latent heating
# above avoids differencing and staggering issues:
# # Ri_moist = similar(Ri); fill!(Ri_e, missing)
# Ri_moist[2:end-1,:] = Ri[2:end-1,:] .* n2moist ./ nm[2:end-1,:]

# Get Python classes for log scale handling
LogNorm = pyimport("matplotlib.colors").LogNorm
LogFormatterMathtext = pyimport("matplotlib.ticker").LogFormatterMathtext

# plots
# timeheight log10 vapor, log10 cloud, theta
clf()
levels = 10 .^ (-6:0.25:-1.5)
cs = contourf(t/24, zs, qv, levels=levels, norm=LogNorm(), cmap=ColorMap("RdYlBu_r"))
#colorbar()
colorbar(cs, format=LogFormatterMathtext())
contour(t/24, zs, qc, levels=[1e-7, 1e-6], linewidths=0.6, colors=["g"], linestyles="-") #, cmap=ColorMap("RdYlBu_r"))
# colorbar()
contour(t/24, zs, th, colors=["k"], linewidths=0.2, levels=280:2:500)
ylim([0,16])
title("Arabian Sea\n specific humidity and potential temperature")
xticks(0:4)
xlabel("time (d)")
ylabel("height (km)")

fmts = ["svg", "pdf", "eps", "png"]
[savefig("ASqvth.$(f)") for f in fmts]

#=
# timeheight of u, Ri
clf()
contourf(t/24, zs, ds["u"][1,1,:,:], levels=-20:1:20, vmin=-15, vmax=15, cmap=ColorMap("RdYlBu_r"))
colorbar()
contour(t/24, zs, th, colors=["k"], linewidths=0.2, levels=280:2:500)
contour(t/24, zw, Ri, levels=[1, 4], colors="c" )
contour(t/24, zw, Ri, levels=[0, 0.25], colors="m" ) # magenta is unstable
# contour(t/24, zw, Ri_moist, levels=[0.25], colors="g", linewidth=0.1 )
# contour(t, zs, ds["rtke"][1,1,:,:], colors="y")
ylim([0, 8])
title("spin up")
title("Arabian Sea\n wind, Ri")
xlabel("time (d)")
xticks(0:4)
ylabel("height (km)")

# moist process only slightly increases instability, esp at 4 km center of a th ML.

[savefig("ASuthRi.$(f)") for f in fmts]
=#
