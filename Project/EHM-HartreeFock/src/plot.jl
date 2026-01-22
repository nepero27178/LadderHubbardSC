#!/usr/bin/julia

# Arguments handler
if length(ARGS) != 1
	println("How to use this program?
Type the following: \$ julia ./plot.jl --mode
Where:
· mode = --scan / --heatmap / --RMPs")
	exit()
else
	UserInput = ARGS
	Mode = UserInput[1][3:end]
end
InMode::String = Mode

# Includer
PROJECT_ROOT = @__DIR__
PROJECT_ROOT *= "/.."   # Up to the effective root
if Mode=="scan"
	include(PROJECT_ROOT * "/src/setup/scan-simulations-setup.jl")
	include(PROJECT_ROOT * "/src/modules/methods-2D-plotting.jl")
if Mode=="heatmap"
	include(PROJECT_ROOT * "/src/setup/heatmap-simulations-setup.jl")
	include(PROJECT_ROOT * "/src/modules/methods-3D-plotting.jl")
elseif Mode=="RMPs" 
	include(PROJECT_ROOT * "/src/setup/heatmap-simulations-setup.jl")
	include(PROJECT_ROOT * "/src/modules/methods-3D-plotting.jl")
	InMode = "heatmap"
else
	@error "Invalid argument. Use: mode = --scan / --heatmap / --RMPs "
end

function main()

	# Read files
	DirPathIn = PROJECT_ROOT * "/simulations/" * InMode *
		"/Setup=$(Setup)/Phase=$(Phase)/Syms=$(Syms...).csv"

	# Create output directory
	DirPathOut = PROJECT_ROOT * "/analysis/Phase=$(Phase)/" * Mode * "/Setup=$(Setup)/"
	mkpath(DirPathOut)

	# Run plot modules
	if Mode=="scan"
		SavePlot2D(
			FilePathIn,
			DirPathOut;
			xVar="V",
			yVar=HFP,
			pVar="δ",
			cs=:winter
		)
		SavePlot2D(
			FilePathIn,
			DirPathOut;
			xVar="δ",
			yVar=HFP,
			pVar="V",
			cs=:winter
		)
	elseif Mode=="heatmap"
		@error "Under construction!"
	elseif Mode=="RMPs"
		@error "Under construction!"
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
