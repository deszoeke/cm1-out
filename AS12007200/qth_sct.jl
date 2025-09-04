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

# histogram
h = fit(Histogram, (th[:], qv[:]), nbins=(30, 30))
wt = pd(h.weights)

LogNorm = PyPlot.matplotlib.colors.LogNorm
SymLogNorm = PyPlot.matplotlib.colors.SymLogNorm

clf()
subplot(2,2,1)
pcolormesh(h.edges[1], 1e3*h.edges[2], pd(h.weights), 
   cmap = ColorMap("RdYlBu_r"),
   norm = SymLogNorm(linthresh=0.1, vmin=minimum(wt), vmax=maximum(wt)) )

plot(th[1:14:end], 1e3*qv[1:14:end], marker=".", linestyle="none", markersize=0.01, color="k")

fmts = ["svg", "pdf", "eps", "png"]
[ savefig("qth_sct.$(f)") for f in fmts ]

