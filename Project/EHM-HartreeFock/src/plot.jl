#!/usr/bin/julia

# Arguments handler
if length(ARGS) != 1
	println("How to use this program?
Type the following: \$ julia ./plot.jl --mode
Where:
· mode = --scan / --heatmap / --RMPs / --record-g")
	exit()
else
	UserInput = ARGS
	Mode = UserInput[1][3:end]
end
InMode::String = Mode

# Includer
PROJECT_ROOT = @__DIR__
PROJECT_ROOT *= "/.."   # Up to the effective root
if in(Mode, ["scan", "heatmap", "record-g"])
	include(PROJECT_ROOT * "/src/setup/" * Mode * "-simulations-setup.jl")
elseif Mode=="RMPs" 
	include(PROJECT_ROOT * "/src/setup/heatmap-simulations-setup.jl")
	InMode = "heatmap"
else
	@error "Invalid argument. Use: mode = --scan / --heatmap / --RMPs " *
		"/ --Record-g"
end
include(PROJECT_ROOT * "/src/setup/graphic-setup.jl")
include(PROJECT_ROOT * "/src/modules/methods-plotting.jl")

function main()

	# Read files
	DirPathIn = PROJECT_ROOT * "/simulations/" *
		InMode * "/Setup=$(Setup)/Phase="
	if Phase=="SU-Singlet"
		FilePathIn = DirPathIn * "SU-Singlet/Syms=$(Syms...).txt"
	elseif Phase=="SU-Triplet"
		FilePathIn = DirPathIn * "SU-Triplet/Syms=$(Syms...).txt"
	else
		FilePathIn = DirPathIn * Phase * ".txt"
	end

	# Create output directory
	DirPathOut = PROJECT_ROOT * "/analysis/Phase=" * Phase * "/" *
		Mode * "/Setup=$(Setup)/"
	mkpath(DirPathOut)

	# Run plot modules
	if Mode=="scan"
		PlotOrderParameter(
			Phase,
			FilePathIn,
			DirPathOut;
			Syms,
			xVar="V",
			pVar="δ",
			cs=:winter,
			RenormalizeHopping,
			Extension="png"
		)
		PlotOrderParameter(
			Phase,
			FilePathIn,
			DirPathOut;
			xVar="δ",
			pVar="V",
			Skip=2,
			cs=:winter,
			RenormalizeHopping,
			Extension="png"
		)
	elseif Mode=="heatmap"
		PlotOrderParameter2D(
			Phase,
			FilePathIn,
			DirPathOut;
#            xVar="U",
#            yVar="V",
			xVar="V",
			yVar="δ",
			cs=:imola,
			Extension="pdf"
		)
	elseif Mode=="RMPs"
		if Phase=="FakeAF"
			@error "It makes no sense to plot RMPs for Phase=$(Phase)!"
		else
			PlotRMPs(
				Phase,
				FilePathIn,
				DirPathOut;
#				xVar="U",
#				yVar="V",
				xVar="V",
				yVar="δ",
				cs=:imola,
				Extension="png"
			)
		end
	elseif Mode=="record-g"
		PlotRecord(
			Phase,
			DirPathIn,
			DirPathOut;
			rVar="g"
		)
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
