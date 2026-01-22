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

@doc raw"""
function GetLabels(
	Phase::String
)::Dict{String,String}

Returns: LaTeX formatted variable labels.
"""
function GetLabels(
	Phase::String
)::Dict{String,String}

	VarLabels::Dict{String,String} = Dict([
		# Variables
		"t" => "t",
		"U" => "U",
		"V" => "V",
		"Lx" => "L_x",
		"δ" => "\\delta",
		"β" => "\\beta",
		# Other
		"ΔT" => "\\Delta T",
		"I" => "\\text{Total steps}",
		"μ" => "\\mu",
		"g0" => "\\g_0",
		"g" => "g",
		"fMFT" => "f_\\mathrm{MFT}"
	])

	if in(Phase, ["AF", "FakeAF"])
		PhaseLabels::Dict{String,String} = Dict([
			# HFPs
			"m" => "m",
			"w0" => "w^{(\\mathbf{0})}",
			"wp" => "w^{(\\bm{\\pi})}",
			# Convergence
			"Qm" => "Q(m)",
			"Qw0" => "Q(w^{(\\mathbf{0})})",
			"Qwp" => "Q(w^{(\\bm{\\pi})})",
			# RMPs
			"RReΔ" => "\\mathrm{Re}\\{\\tilde{\\Delta}_\\mathbf{k}\\}",
			"RImΔ" => "\\mathrm{Im}\\{\\tilde{\\Delta}_\\mathbf{k}\\}",
			"Rt" => "\\tilde{t}"
		])
	elseif in(Phase, ["SU-Singlet", "FakeSU-Singlet"])
		PhaseLabels = Dict([
			# HFPs
			"Δs" => "|\\Delta^{(s)}|",
			"ΔS" => "|\\Delta^{(s*)}|",
			"Δd" => "|\\Delta^{(d)}|",
			"gS" => "g^{(s*)}",
			"gd" => "g^{(d)}",
			# Convergence
			"QΔs" => "Q(\\Delta^{(s)})",
			"QΔS" => "Q(\\Delta^{(s*)})",
			"QΔd" => "Q(\\Delta^{(d)})",
			"QgS" => "Q(g^{(s*)})",
			"Qgd" => "Q(g^{(d)})",
			# RMPs
			"t_tilde" => "\\tilde{t}"
		])
	elseif in(Phase, ["SU-Triplet", "FakeSU-Triplet"])
		PhaseLabels = Dict([
			# HFPs
			"Δpx" => "|\\Delta_{p_x}",
			"Δpy" => "\\Delta_{p_y}",
			"Δp+" => "\\Delta_{p_+}",
			"Δp-" => "\\Delta_{p_-}",
			"gS" => "g^{(s^*)}",
			"gd" => "g^{(d)}",
			# Convergence
			"QΔpx" => "Q(|\\Delta_{p_x})",
			"QΔpy" => "Q(\\Delta_{p_y})",
			"QΔp+" => "Q(\\Delta_{p_+})",
			"QΔp-" => "Q(\\Delta_{p_-})",
			"QgS" => "Q(g^{(s^*)})",
			"Qgd" => "Q(g^{(d)})",
			# RMPs
			"t_tilde" => "\\tilde{t}"
		])
	end

	return merge(VarLabels, PhaseLabels)

end

function Plot2D(
	FilePathIn::String;					# Data filepath
	Print::Bool=false,					# Set true to format for savefig
	xVar::String="V",					# Specify x variable
	yVar::String="m",					# Specify y variable
	pVar::String="δ",					# Specify parametric variable
	cs::Symbol=:imola50,					# Custom colorscheme
	Skip::Int64=0,						# pVar skip parameter
)::Vector{GroupedPlot}

	# List vars and pars
	xVars = ["t", "U", "V", "Lx", "δ", "β"]
	Pars = filter(!=(pVar), filter(!=(xVar), xVars))

	# Input safecheck
	!in(xVar, xVars) ? error("Invalid x variable, choose one of $(xVars)") : 0
	!in(pVar, xVars) ? error("Invalid p variable, choose one of $(xVars)") : 0
	xVar==pVar ? error("You have chosen xVar=pVar!") : 0

	# Unpack filepath
	Setup, Phase, Syms = UnpackFilePath(FilePathIn)
	RenormalizeBands::Bool=true
	occursin("Fake",Phase) ? RenormalizeBands=false : 0

	# Load data
	DF = ImportData(FilePathIn)
	Sim::Simulation = Simulation(DF,Setup,Phase,Syms)
	DF = EnlargeDF!(Sim) # Unpack dictionaries, compute RMPs
	yVars = filter(!in(xVars), names(DF))
	!in(yVar,yVars) ? error("Invalid y variable, choose one of $(yVars)") : 0

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
	GroupedDF = groupby(DF,Pars) #TODO Add skip
	J = length(GroupedDF)
	C = floor(Int64, length(colorschemes[cs]) / J)
	PlotVec = GroupedPlot[]

	# Cycle over simulated points
	for (j,df) in enumerate(GroupedDF)

		# Select plot parameters and print on terminal
		PltPars = DataFrame(select(df, Symbol.(Pars))[1,:])
		InfoPars = copy(PltPars)
		InfoPars[!,"x"] .= xVar
		InfoPars[!,"y"] .= yVar
		InfoPars[!,"p"] .= pVar
		@info "\e[1;36mScan plot $(j)/$(J)\e[0m" InfoPars
		println()

		yLabel::String = ""
		if yVar=="m"
			yLabel = "Magnetization"
		elseif Print
			yLabel = "\$" * VarLabels[yVar] * "\$"
		elseif !Print
			yLabel = yVar
		end

		# Initialize plot
		Fig = Figure(size=(600,400),figure_padding = 1)
		ax = Axis(Fig[1, 1])
		if Print
			ax.xlabel = L"$%$(VarLabels[xVar])$"
			ax.ylabel = L"$%$(VarLabels[yVar])$"
		elseif !Print
			ax.xlabel = xVar
			ax.ylabel = yVar
		end

		# Create filename
		FileName = join( ["$(Pars[i])=$(df[!,Pars[i]][1])" for i in 1:length(Pars)], '_' )
		FileName = yVar * "_" * FileName

		# Create raw title string, either for printing or local plotting
		rawTitle::String = yLabel * " ("
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

		Groupeddf = groupby(df,pVar)
		I = length(Groupeddf)
		C = floor(Int64, length(colorschemes[cs]) / I)
		C==0.0 ? error("Your cs (ColorScheme) is not large enough.") : 0
		for (i,pdf) in enumerate(Groupeddf)
			p = pdf[!,pVar][1]
			xx = pdf[!,xVar]
			yy = pdf[!,yVar]

			if Print
				label = L"$%$(VarLabels[pVar])=%$(p)$"
			elseif !Print
				label = pVar * "=$(p)"
			end

			# Plot parametrically
			in(yVar, GetHFPs(Phase;Syms)) ? yy = abs.(yy) : 0
			scatter!(
				ax, xx, yy,
				marker = :circle,
				color = colorschemes[cs][C*i],
				markersize = 8,
				label = label
			)
			lines!(
				ax, xx, yy,
				color = colorschemes[cs][C*i]
			)
		end

		if Print
			LegendLabel = L"Simulated $%$(pVar)$:"
		elseif !Print
			LegendLabel = "Simulated " * pVar * ":"
		end

		MinX = minimum(df[!,xVar])
		MaxX = maximum(df[!,xVar])
		MinY = minimum(df[!,yVar])
		MaxY = maximum(df[!,yVar])

		Pad = 5e-2
		MinXLim = MinX - Pad*(MaxX-MinX)
		MaxXLim = MaxX + Pad*(MaxX-MinX)
		MinYLim = MinY - Pad*(MaxY-MinY)
		MaxYLim = MaxY + Pad*(MaxY-MinY)

		xlims!(ax, MinXLim, MaxXLim)
		ylims!(ax, MinYLim, MaxYLim)
		Fig[1, 2] = Legend(Fig, ax, LegendLabel, framevisible = false)
		push!(PlotVec, GroupedPlot(Fig,df,FileName))
	end

	return PlotVec
end

function SavePlot2D(
	FilePathIn::String,					# Data filepath
	DirPathOut::String;					# Output directory
	xVar::String="V",					# Specify x variable
	yVar::String="m",					# Specify y variable
	pVar::String="δ",					# Specify parametric variable
	cs::Symbol=:imola50,					# Custom colorscheme
	Extension::String="pdf"				# File extension
)

	# Assert printing
	Print::Bool=true
	PlotVec = Plot2D(FilePathIn;xVar,yVar,pVar,Print)

	# Initialize directory structure
	Setup, Phase, Syms = UnpackFilePath(FilePathIn)
	DirPathOut *= "/Syms=$(Syms...)/xVar=" * xVar * "_pVar=" * pVar * "/"
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