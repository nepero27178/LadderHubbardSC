#!/usr/bin/julia
PROJECT_ROOT = @__DIR__
PROJECT_ROOT *= "/.."
include(PROJECT_ROOT * "/src/setup/simulations-setup.jl")
include(PROJECT_ROOT * "/src/setup/graphic-setup.jl")
include(PROJECT_ROOT * "/src/modules/methods-plotting.jl")

function main()
    DirPathIn = PROJECT_ROOT * "/simulations/hubbard/p=$(p)_Δm=$(Δm)_Δn=$(Δn)/"
    FilePathIn = DirPathIn * "/Setup=$(Setup).txt"
    DirPathOut = PROJECT_ROOT * "/analysis/hubbard/p=$(p)_Δm=$(Δm)_Δn=$(Δn)/"
    mkpath(DirPathOut)
    PlotUΔ(FilePathIn, DirPathOut)  # Δ vs U
    PlotδΔ(FilePathIn, DirPathOut)  # Δ vs δ
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
