#!/usr/bin/julia

using GLMakie
using CairoMakie
using LaTeXStrings
using ColorSchemes
using DataFrames
using DelimitedFiles

PROJECT_METHODS_DIR = @__DIR__
include(PROJECT_METHODS_DIR * "/methods-IO.jl")
include(PROJECT_METHODS_DIR * "/structs.jl")

function Plot3D(
	FilePathIn::String;					# Data filepath
	Print::Bool=false,					# Set true to format for savefig
	Mode::String="heatmap",				# Specify plot method
	xVar::String="U",					# Specify x variable
	yVar::String="V",					# Specify y variable
	zVar::String="fMFT",					# Specify z variable
	cs::Symbol=:imola50,					# Custom colorscheme
)::Vector{GroupedPlot}

	# Define plot function (heatmap, surface, mesh)
	PlotFunction = eval(Meta.parse("CairoMakie." * Mode * "!"))

	# List vars and pars
	xyVars = ["t", "U", "V", "Lx", "δ", "β"]
	Pars = filter(!=(yVar), filter(!=(xVar), xyVars))

	# Input safecheck
	!in(xVar, xyVars) ? error("Invalid x variable, choose one of $(xyVars)") : 0
	!in(yVar, xyVars) ? error("Invalid y variable, choose one of $(xyVars)") : 0
	xVar==yVar ? error("You have chosen xVar=yVar!") : 0

	# Unpack filepath
	Setup, Phase, Syms = UnpackFilePath(FilePathIn)
	RenormalizeBands::Bool=true
	occursin("Fake",Phase) ? RenormalizeBands=false : 0

	# Load data
	DF = ImportData(FilePathIn)
	Sim::Simulation = Simulation(DF,Setup,Phase,Syms)
	EnlargeDF!(Sim) # Unpack dictionaries, compute RMPs
	zVars = filter(!in(xyVars), names(DF))
	!in(zVar,zVars) ? error("Invalid z variable, choose one of $(zVars)") : 0

	if Print
		# Activate backend
		CairoMakie.activate!()

		MT = Makie.MathTeXEngine
		MT_DIR = dirname(pathof(MT)) * "/../assets/fonts/NewComputerModern"

		set_theme!(fonts = (
			regular = MT_DIR * "/NewCM10-Regular.otf",
			bold = MT_DIR * "/NewCM10-Bold.otf"
		))

		# Get LaTeX formatted labels
		VarLabels = GetLabels(Phase)
	elseif !Print
		# Activate backend
		GLMakie.activate!()
	end

	# Group data
	GroupedDF = groupby(DF,Pars)
	J = length(GroupedDF)
	PlotVec = GroupedPlot[]

	# Cycle over simulated points
	for (j,df) in enumerate(GroupedDF)

		# Select plot parameters and print on terminal
		PltPars = DataFrame(select(df, Symbol.(Pars))[1,:])
		InfoPars = copy(PltPars)
		InfoPars[!,"x"] .= xVar
		InfoPars[!,"y"] .= yVar
		InfoPars[!,"z"] .= zVar
		@info "\e[1;36mHeatmap plot $(j)/$(J)\e[0m" InfoPars
		println()

		# Reshape local data
		xx, yy, zz = ReshapeData(DataFrame(df);xVar,yVar,zVar)

		zLabel::String = ""
		if zVar=="m"
			zLabel = "Magnetization"
		elseif Print
			zLabel = "\$" * VarLabels[zVar] * "\$"
		elseif !Print
			zLabel = zVar
		end

		# Initialize plot
		H = Figure(size=(600,400),figure_padding = 1)
		ax = Axis(H[1, 1])
		if Print
			ax.xlabel = L"$%$(VarLabels[xVar])$"
			ax.ylabel = L"$%$(VarLabels[yVar])$"
		elseif !Print
			ax.xlabel = xVar
			ax.ylabel = yVar
		end

		# Create filename
		FileName = join( ["$(Pars[i])=$(df[!,Pars[i]][1])" for i in 1:length(Pars)], '_' )
		FileName = zVar * "_" * FileName

		# Create raw title string, either for printing or local plotting
		rawTitle::String = zLabel * " ("
		if Print
			ParTitle = [VarLabels[Par] * "=$(PltPars[!,Par][1])" for Par in Pars]
			ParTitle = ["\$" * Par * "\$" for Par in ParTitle]
		elseif !Print
			ParTitle = [Par * "=$(PltPars[!,Par][1])" for Par in Pars]
		end
		rawTitle *= join(ParTitle, ", ") * ")"

		# Include RenormalizeBands specifications
		if !RenormalizeBands && Print
			r = split(rawTitle, "t=")
			rawTitle = r[1] * "t=\\tilde{t}=" * r[2]
		elseif !Print
			r = split(rawTitle, ")")
			rawTitle = r[1] * ", rb=$(RenormalizeBands))"
		end

		# Handle infinities
		Print ? rawTitle = replace(rawTitle, "Inf" => "\\infty") : 0
		if Print
			ax.title = L"%$(rawTitle)"
		elseif !Print
			ax.title = rawTitle
		end

		# Plot parametrically
		in(zVar, GetHFPs(Phase;Syms)) ? zz = abs.(zz) : 0
		h = PlotFunction(
			xx, yy, zz',
			colormap=cs
		)
		Colorbar(H[1,2], h)

		push!(PlotVec, GroupedPlot(H,df,FileName))
	end

	return PlotVec
end

function SavePlot3D(
	FilePathIn::String,					# Data filepath
	DirPathOut::String;					# Output directory
	Mode::String="heatmap",				# Specify plot method
	xVar::String="U",					# Specify x variable
	yVar::String="V",					# Specify y variable
	zVar::String="fMFT",					# Specify z variable
	cs::Symbol=:imola50,					# Custom colorscheme
	Extension::String="pdf"				# File extension
)

	# Assert printing
	Print::Bool=true
	PlotVec = Plot3D(FilePathIn;Mode,xVar,yVar,zVar,Print)

	# Initialize directory structure
	Setup, Phase, Syms = UnpackFilePath(FilePathIn)
	DirPathOut *= "/xVar=" * xVar * "_yVar=" * yVar * "/"
	mkpath(DirPathOut)

	# Save each plot
	for GP in PlotVec
		FilePathOut = DirPathOut * "/" * GP.FileName * "." * Extension
		with_theme(theme_latexfonts()) do #TODO Redundant?
			save(FilePathOut,GP.H.scene)
		end
	end

	@info "\e[1;36mPlots saved at:\e[0m" DirPathOut

end