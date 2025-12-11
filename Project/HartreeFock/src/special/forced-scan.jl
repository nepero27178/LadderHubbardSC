#!/usr/bin/julia
using DelimitedFiles
using DataFrames

# Includer
PROJECT_ROOT = @__DIR__
PROJECT_ROOT *= "/../.." # Up to the effective root
include(PROJECT_ROOT * "/src/setup/graphic-setup.jl")
include(PROJECT_ROOT * "/src/modules/methods-simulating.jl")
include(PROJECT_ROOT * "/src/modules/methods-plotting.jl")

function RunHFScan_AF_ForceHoppingShift(
	DF::DataFrame,						# Input DataFrame to be reproduced
	FilePathOut::String="",				# Output file
)

	Phase = "AF"
	p = 100
	Δv = Dict([
		"m" => 1e-4,
		"w0" => 1e-4,
		"wp" => 1e-4
	])
	Δn = 1e-2
	g = 0.5

	# Get Hartree Fock Parameters labels
    HFPs = GetHFPs(Phase)

	# File coditional initialization (otherwise, just append)
	if FilePathOut != ""
		Header = "t;U;V;Lx;β;δ;v;Q;ΔT;μ\n"
		write(FilePathOut, Header)
	end

	# HF iterations
	i = 1
	Iterations = size(DF,1)
	for Point in eachrow(DF)

		# Initializer
		v0::Dict{String,Float64} = eval(Meta.parse(Point.v))

		# Model parameters
		w0 = eval(Meta.parse(Point.v))["w0"]
		t = Point.t - w0 * Point.V
		U = Point.U
		V = Point.V

		Parameters::Dict{String,Float64} = Dict([
			"t" => t,
			"U" => U,
			"V" => V
		])

		Lx = Int64(Point.Lx)
		β = Point.β
		δ = Point.δ
		L = [Lx, Lx]
		printstyled(
			"\e[2K\e[1GRun ($i/$Iterations): " *
			"$Phase HF at t=$t, U=$U, V=$V, L=$Lx, β=$β, δ=$δ",
			color=:yellow
		)
		ResultsVector::Matrix{Any} = [0 0]  # Dummy placeholder
		ResultsVector = hcat(ResultsVector, [t U V Lx β δ])

		# Run routine, all positional arguments here must be false
		HFResults = RunHFAlgorithm(
			Phase,Parameters,L,0.5+δ,β,
			p,Δv,Δn,g;
			v0i=v0,
			RenormalizeHopping=false
		)

		v::Dict{String,Float64} = Dict([
			key => HFResults[1][key] for key in HFPs
		])
		Qs::Dict{String,Float64} = Dict([
			key => HFResults[2][key] for key in HFPs
		])
		ResultsVector = hcat(ResultsVector[:,3:end], [v Qs HFResults[3] HFResults[4]])

		i += 1

		# Append to initialized or existing file
		if FilePathOut != ""
			open(FilePathOut, "a") do io
				writedlm(io, ResultsVector, ';')
			end
		end

		v0 = copy(v)
	end

	printstyled(
		"\e[2K\e[1GDone! Data saved at " * FilePathOut * "\n", color=:green
	)

end

function ForcedScan(
	Setup::String,
)
	FilePathIn = PROJECT_ROOT * "/simulations/scan/Setup=$(Setup)/AF.txt"
	DataIn = open(FilePathIn) do io
		readdlm(FilePathIn, ';', comments=false, '\n')
	end
	DF = DataFrame(DataIn[2:end,:], DataIn[1,:])

	FilePathOut = replace(FilePathIn, "scan" => "special/forced-scan")
	mkpath(dirname(FilePathOut))

	RunHFScan_AF_ForceHoppingShift(DF,FilePathOut)
end

function PlotForcedScan(
	Setup::String;
	xVar::String="U",					# Specify x variable
	pVar::String="V",					# Specify parametric variable
	Skip::Int64=1,						# xVar skip parameter
	cs::Symbol=:imola50,				# Custom colorscheme
	Extension::String="pdf"				# File extension
)

	FilePathIn = PROJECT_ROOT * "/simulations/special/forced-scan/Setup=$(Setup)/AF.txt"
	DirPathOut = PROJECT_ROOT * "/analysis/Phase=AF/special/forced-scan/Setup=$(Setup)/"

	Phase = "AF"
	RenormalizeHopping = true
	FilePathOut = ""
	AllVars = ["t", "U", "V", "Lx", "δ", "β"]

	# Input safecheck
	if !in(xVar, AllVars)
		@error "Invalid x variable, choose one of $(AllVars)"
		exit()
	end

	if !in(pVar, AllVars)
		@error "Invalid parametric variable, choose one of $(AllVars)"
		exit()
	end

	if xVar==pVar
		@error "You have chosen xVar=pVar!"
		exit()
	end

	# Get LaTeX formatted labels
	VarLabels = GetLabels(Phase)

	# Initialize directory structure
	DirPathOut *= "PlotOrderParameter/xVar=" * xVar * "/pVar=" * pVar * "/"
	mkpath(DirPathOut)

	# Read data and create DataFrame
	DataIn = open(FilePathIn) do io
		readdlm(FilePathIn, ';', comments=false, '\n')
	end
	DF = DataFrame(DataIn[2:end,:], DataIn[1,:])
	DF.Lx = Int64.(DF.Lx)

	# Unique variables dictionary
	uDict::Dict{String,Vector{Any}} = Dict([
		Var => unique(DF[!, Symbol(Var)]) for Var in AllVars
	])

	uDict[xVar] = [NaN]					# Remove from following cycle
	JJ::Vector = uDict[pVar][1:Skip:end]
	uDict[pVar] = [NaN]					# Remove from following cycle
	q = floor(Int64, length(colorschemes[cs]) / length(JJ) )

	# Check if the selected colorscheme is large enough
	if q==0
		@error "Your cs (ColorScheme) is not large enough. Select another."
		exit()
	end

	# List HF parameters
	ListHFPs::Vector{String} = [key for key in keys(
		eval(Meta.parse(
			DF.v[1]
		))
	)]

	# Cycle over HF parameters
	for HF in ListHFPs

		HFLabel::String = "\$" * VarLabels[HF] * "\$"
		if HF=="m"
			HFLabel = "Magnetization"
		end

		# Cycle over simulated points
		for (u,U) in enumerate(uDict["U"]),
			(v,V) in enumerate(uDict["V"]),
			(l,L) in enumerate(uDict["Lx"]),
			(d,δ) in enumerate(uDict["δ"]),
			(b,β) in enumerate(uDict["β"])

			# Initialize local dictionary
			lDict::Dict{String,Any} = Dict([
				#"t" => t,
				"U" => U,
				"V" => V,
				"Lx" => L,
				"δ" => δ,
				"β" => β
			])

			# Select entries
			Selections::Vector{Bool} = [true for _ in 1:length(DF.t)]
			for Var in AllVars
				if !in(Var, ["t", xVar, pVar])
					"""
					Selections = Selections .* (DF[Var] .== lDict[Var])
					"""
					Selections = Selections .* (DF[!, Symbol(Var)] .== lDict[Var])
				end
			end

			# Initialize plot
			P = plot(
				size = (600,400),
				xlabel = L"$%$(VarLabels[xVar])$",
				ylabel = L"$%$(VarLabels[HF])$",
				legend = :outertopright,
			)

			# Write terminal message and file name
			TerminalMsg::String = "\e[2K\e[1GPlotting $(Phase) HF $(HF) data: "
			FilePathOut::String = DirPathOut * "/" * HF
			rawTitle::String = HFLabel * " ("

			AvoidList::Vector{String} = [xVar, pVar] #TODO Extend to pure hubbard
			for Var in AllVars
				if !in(Var, AvoidList)
					if Var=="t"
						rawTitle *= "\$t=1.0-w^{(\\mathbf{0})}V\$, "
					else
						lVar = lDict[Var]
						TerminalMsg *= Var * "=$(lVar), "
						FilePathOut *= "_" * Var * "=$(lVar)"
						rawTitle *= "\$$(VarLabels[Var])=$(lVar)\$, "
					end
					
					if Var=="Lx" # I am desperate about correct formatting
						rawTitle = rawTitle[1:end-5] * "\$, "
					elseif Var=="β" && β==Inf
						rawTitle = rawTitle[1:end-6] * "\\infty\$, "
					end
				end
			end
			TerminalMsg = TerminalMsg[1:end-2] * " [x variable: " * xVar *
				", parametric variable: " * pVar * "]"
			FilePathOut *= "." * Extension
			rawTitle *= "varying \$" * VarLabels[pVar] * "\$)"
			printstyled("\e[2K\e[1G" * TerminalMsg, color=:yellow)

			# Cycle over parametric variable
			for (j,J) in enumerate(JJ)

				# Run local selection
				lSelections = Selections .* (DF[!, Symbol(pVar)] .== J)
				vv = DF.v[lSelections]

				# Define x and y variables
				xx::Vector{Float64} = DF[!, Symbol(xVar)][lSelections]
				yy::Vector = [eval(Meta.parse(
					vv[i]
				))[HF] for i in 1:length(vv)]

				# Plot parametrically
				plot!(
					xx, yy,
					markershape = :circle,
					markercolor = colorschemes[cs][q*j],
					markersize = 1.5,
					linecolor = colorschemes[cs][q*j],
					label = L"$%$(VarLabels[pVar])=%$(J)$",
					legendfonthalign = :left
				)
				title!(L"%$(rawTitle)")
			end

			# Save figure
			Plots.PGFPlotsX.push_preamble!(
				backend_object(P).the_plot, "\\usepackage{bm}"
			)
			savefig(P, FilePathOut)
		end
	end
	printstyled("\e[2K\e[1GDone! Plots saved at $(DirPathOut)\n", color=:green)
end

if abspath(PROGRAM_FILE) == @__FILE__
	Setup = "B[256]"
	#ForcedScan(Setup)
	PlotForcedScan(
		Setup;
		xVar="V", 
		pVar="δ",
		cs=:winter,
		Extension="png"
	)
end