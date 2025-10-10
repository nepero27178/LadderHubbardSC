#!/usr/bin/julia
PROJECT_ROOT = @__DIR__
PROJECT_ROOT *= "/.."
include(PROJECT_ROOT * "/src/setup/simulations-setup.jl")
include(PROJECT_ROOT * "/src/setup/graphic-setup.jl")
include(PROJECT_ROOT * "/src/modules/methods-plotting.jl")

function main()
    DirPathIn = PROJECT_ROOT * "/simulations/Setup=$(Setup)/"
    FilePathIn = DirPathIn * "$(HMSymsStr)-wave.txt"
    DirPathOut = PROJECT_ROOT * "/analysis/"
    mkpath(DirPathOut)
    # PlotVΔ(HMSyms, FilePathIn, DirPathOut)  # Δ vs V
    # PlotδΔ(HMSyms, FilePathIn, DirPathOut)  # Δ vs δ
	PlotHeatmapVdΔ(HMSyms, FilePathIn, DirPathOut)
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
