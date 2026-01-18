using Plots
using DataFrames

PROJECT_LABS_DIR = @__DIR__
include(PROJECT_LABS_DIR * "/../../modules/methods-optimizations.jl")

L = [50, 50]
Kx::Vector{Float64} = [kx for kx in -1:2/L[1]:1]
popfirst!(Kx)
Ky::Vector{Float64} = [ky for ky in -1:2/L[2]:1]
popfirst!(Ky)
K::Matrix{Vector{Float64}} = [ [kx,ky] for kx in Kx, ky in Ky ]

Sym = "S"
wk = GetWeight.(K;Sym)
DF = DataFrame(
    x = [k[1] for k in K][:],
    y = [k[2] for k in K][:],
    z = wk[:]
)

DF0 = filter(row -> row.z == 0, DF)
DF1 = filter(row -> row.z == 1, DF)
DF2 = filter(row -> row.z == 2, DF)
DF4 = filter(row -> row.z == 4, DF)

scatter(
    DF0.x, DF0.y,
    color="gray75",
    markersize=1,
    markeralpha=0,
    markerstrokewidth=0,
    label="w=0",
    legend=:outertopleft
)
scatter!(
    DF1.x, DF1.y,
    color="green",
    markersize=2,
    markerstrokewidth=0.5,
    label="w=1"
)
scatter!(
    DF2.x, DF2.y,
    color="yellow",
    markersize=2,
    markerstrokewidth=0.5,
    label="w=2"
)
scatter!(
    DF4.x, DF4.y,
    color="cyan",
    markersize=2,
    markerstrokewidth=0.5,
    label="w=4"
)
xlabel!("kx/pi")
ylabel!("kx/pi")
title!("Weights for sym: $(Sym)")
