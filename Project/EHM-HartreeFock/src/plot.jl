#!/usr/bin/julia

using GLMakie
using CairoMakie
using LaTeXStrings
using ColorSchemes
using DataFrames
using DelimitedFiles

# Includer
PROJECT_SRC_DIR = @__DIR__
include(PROJECT_SRC_DIR * "/modules/methods-IO.jl")

# Arguments handler
if length(ARGS) != 2
	println("How to use this program?
Type the following: \$ julia plot.jl --mode --obj
Where:
· mode = hs / rs
· obj = HFPs / RMPs / Qs / phys")
	exit()
else
	UserInput = ARGS
	Mode = UserInput[1][3:end]
	Obj = UserInput[2][3:end]
end

# Process mode
if in(Mode, ["rs","hs"])
	include(PROJECT_SRC_DIR * "/setup/" * Mode * "-setup.jl")
	Mode=="hs" ? Dim = "2D" : Dim = "3D"
	include(PROJECT_SRC_DIR * "/modules/methods-" * Dim * "-plotting.jl")
else
	@error "Invalid mode. Use: mode = hs / rs"
end

# Process obj
if Obj=="HFPs"
	objList = GetHFPs(Phase;Syms)
elseif Obj=="RMPs"
	objList = [key for key in keys(GetRMPs(Phase;Syms))]
elseif Obj=="Qs"
	QsList = ["Q$(HFP)" for HFP in GetHFPs(Phase;Syms)]
	RunList = ["ΔT", "I", "g0", "g"]
	objList = vcat(QsList, RunList)
elseif Obj=="phys"
	objList = ["μ", "fMFT"]
else
	@error "Invalid obj. Use obj = HFPs / RMPs / Qs / phys"
end

function main()

	# Read files
	FilePathIn = PROJECT_SRC_DIR * "/../simulations/Mode=$(Mode)/Setup=$(Setup)/Phase=$(Phase)/Syms=$(Syms...).csv"

	# Create output directory
	# For simulations: Setup > Phase > Syms (to make comparable data in the same folder)
	# For plots: Phase > Setup > Syms (to make same-phase plots in the same folder)
	DirPathOut = PROJECT_SRC_DIR * "/../analysis/Mode=$(Mode)/Phase=$(Phase)/Setup=$(Setup)/Syms=$(Syms...)/Obj=$(Obj)"
	mkpath(DirPathOut)

	for obj in objList
		# Run plot modules for each HFP
		if Mode=="hs"
			SavePlot2D(
				FilePathIn,
				DirPathOut;
				xVar="V",
				yVar=obj,
				pVar="δ",
				# cs=:winter
			)
			SavePlot2D(
				FilePathIn,
				DirPathOut;
				xVar="δ",
				yVar=obj,
				pVar="V",
				# cs=:winter
			)
		elseif Mode=="rs"
			SavePlot3D(
				FilePathIn,
				DirPathOut;
				xVar="U",
				yVar="V",
				zVar=obj,
				# cs=:winter
			)
		end
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
