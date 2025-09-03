# plot CM1 output

activate("..")

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

# load surface slice from 3D data
ds = NCDataset("data/cm1out_s.nc")
lat = ds["lat"][:]
lon = ds["lon"][:]
lev = ds["lev"][:]
time = ds["time"][:]
gt(v) = ds[v][:,:,1,5]
th = gt("th")
qv = gt("qv")
u10 = ds["u10"][:,:,5]
v10 = ds["v10"][:,:,5]

clf()
subplot(2,1,1)
pcolormesh((0:length(lon)).*0.5, (0:length(lat)).*0.5, pd(th), cmap=ColorMap("RdYlBu_r"))
colorbar(label=L"\theta"*" (K)")
gca().set_aspect("equal")

subplot(2,1,2)
pcolormesh((0:length(lon)).*0.5, (0:length(lat)).*0.5, 1e3.*pd(qv), cmap=ColorMap("RdYlBu_r"))
colorbar(label=L"q_v"*" (g/kg)")
gca().set_aspect("equal")
tight_layout()

fmts = ["svg", "pdf", "eps", "png"]
[ savefig("coldpools_d4.$(f)") for f in fmts ]

clf()
subplot(2,1,1)
pcolormesh((0:length(lon)).*0.5, (0:length(lat)).*0.5, pd(u10), vmin=0, vmax=14, cmap=ColorMap("RdYlBu_r"))
colorbar(label=L"u_{10}"*" (m/s)")
gca().set_aspect("equal")

subplot(2,1,2)
pcolormesh((0:length(lon)).*0.5, (0:length(lat)).*0.5, pd(v10), cmap=ColorMap("RdYlBu_r"))
colorbar(label=L"v_{10}"*" (m/s)")
gca().set_aspect("equal")
tight_layout()

[savefig("u10v10.$(f)") for f in fmts]
