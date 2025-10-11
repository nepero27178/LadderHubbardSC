#!/usr/bin/julia
PROJECT_ROOT = @__DIR__
PROJECT_ROOT *= "/.."
include(PROJECT_ROOT * "/src/setup/simulations-setup.jl")
include(PROJECT_ROOT * "/src/setup/graphic-setup.jl")
include(PROJECT_ROOT * "/src/modules/methods-plotting.jl")

function main()
    DirPathIn = PROJECT_ROOT * "/simulations/Phase=$(Phase)/Setup=$(Setup)/"
    FilePathIn = DirPathIn * "$(Model).txt"
    DirPathOut = PROJECT_ROOT * "/analysis/Phase=$(Phase)/Setup=$(Setup)/"
    mkpath(DirPathOut)
    PlotOrderParameter(
        Phase,
        FilePathIn,
        DirPathOut;
        pVar="V",
        xVar="δ"
    )
	# PlotHeatmapVdΔ(HMSyms, FilePathIn, DirPathOut)
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
