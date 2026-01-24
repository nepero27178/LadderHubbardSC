@doc raw"""
function PlotOrderParameter(
	Phase::String,
	FilePathIn::String,
	DirPathOut::String;
	xVar::String=\"U\",
	pVar::String=\"V\",
	Skip::Int64=1,
	cs::Symbol=:imola50,
	RenormalizeBands::Bool=true,
	Extension::String="pdf"::String="pdf"
)

Returns: none (plots saved at `DirPathOut`).

`PlotOrderParameter` takes as input `Phase` (string specifying the mean-field
phase, the allowed are \"AF\", \"FakeAF\", \"SU-Singlet\", \"SU-Triplet\"),
`FilePathIn` (path to the data files), `DirPathOut` (path to the output
directory). The optional parameters are `xVar` and `pVar` (strings specifying
respectively the x variable and the parametric variable of the plot, the
allowed are \"t\", \"U\", \"V\", \"Lx\", \"δ\", \"β\"), `Skip` (integer
indicating how many steps in the parametric variable to be skipped between two
plots), `cs` (colorscheme symbol). The boolean option `RenormalizeBands'
allows for choosing to renormalize or not the hopping parameter. `Extension`
selects the file extension (the allowed are \"pdf\", \"png\", \"svg\").
"""
function PlotOrderParameter(
	Phase::String,						# Mean field phase
	FilePathIn::String,					# Data filepath
	DirPathOut::String;					# Output directory path
	Syms::Vector{String}=["s"],			# Gap function symmetries
	xVar::String="U",					# Specify x variable
	pVar::String="V",					# Specify parametric variable
	Skip::Int64=1,						# xVar skip parameter
	cs::Symbol=:imola50,					# Custom colorscheme
	RenormalizeBands::Bool=true,		# Conditional renormalization of t
	Extension::String="pdf"				# File extension
)

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
	DirPathOut *= "PlotOrderParameter"
	if in(Phase, ["SU-Singlet", "FakeSU-Singlet"])
		DirPathOut *= "_Syms=$(Syms...)"
	elseif in(Phase, ["SU-Triplet", "FakeSU-Triplet"])
		DirPathOut *= "_Syms=$(Syms...)"
	end
	DirPathOut *= "/xVar=" * xVar * "/pVar=" * pVar * "/"
	mkpath(DirPathOut)

	# Import data
	DF = ImportData(FilePathIn)

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
		for (w,t) in enumerate(uDict["t"]),
			(u,U) in enumerate(uDict["U"]),
			(v,V) in enumerate(uDict["V"]),
			(l,L) in enumerate(uDict["Lx"]),
			(d,δ) in enumerate(uDict["δ"]),
			(b,β) in enumerate(uDict["β"])

			# Initialize local dictionary
			lDict::Dict{String,Any} = Dict([
				"t" => t,
				"U" => U,
				"V" => V,
				"Lx" => L,
				"δ" => δ,
				"β" => β
			])

			# Select entries
			Selections::Vector{Bool} = [true for _ in 1:length(DF.t)]
			for Var in AllVars
				if !in(Var, [xVar, pVar])
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
					lVar = lDict[Var]
					TerminalMsg *= Var * "=$(lVar), "
					FilePathOut *= "_" * Var * "=$(lVar)"
					RS::String=""
					if !RenormalizeBands && Var=="t"
						RS *= "\\tilde{t}="
					end
					rawTitle *= "\$$(VarLabels[Var])=" * RS * "$(lVar)\$, "
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

@doc raw"""
function PlotOrderParameter2D(
	Phase::String,
	FilePathIn::String,
	DirPathOut::String;
	xVar::String=\"U\",
	yVar::String=\"V\",
	cs::Symbol=:imola50,
	Extension::String="pdf"
)

Returns: none (plots saved at `DirPathOut`).

`PlotOrderParameter2D` takes as input `Phase` (string specifying the mean-field
phase, the allowed are \"AF\", \"FakeAF\", \"SU-Singlet\", \"SU-Triplet\"), `FilePathIn`
(path to the data files), `DirPathOut` (path to the output directory). The
optional parameters are `xVar` and `yVar` (strings specifying respectively the x
variable and the y variable of the plot, the allowed are \"t\", \"U\", \"V\",
\"Lx\", \"δ\", \"β\"), `cs` (colorscheme symbol). `Extension` selects the file
extension (the allowed are \"pdf\", \"png\", \"svg\").
"""
function PlotOrderParameter2D(
	Phase::String,						# Mean field phase
	FilePathIn::String,					# Data filepath
	DirPathOut::String;					# Output directory path
	xVar::String="U",                   # Specify x variable
	yVar::String="V",                   # Specify y variable
	cs::Symbol=:imola50,                # Custom colorscheme
	Extension::String="pdf"             # File extension
)

	FilePathOut = ""
	AllVars = ["t", "U", "V", "Lx", "δ", "β"]

	# Input safecheck
	if !in(xVar, AllVars)
		@error "Invalid x variable, choose one of $(AllVars)"
		exit()
	end

	if !in(yVar, AllVars)
		@error "Invalid y variable, choose one of $(AllVars)"
		exit()
	end

	if xVar==yVar
		@error "You have chosen xVar=yVar!"
		exit()
	end

	# Get LaTeX formatted labels
	VarLabels = GetLabels(Phase)

	# Initialize directory structure
	DirPathOut *= "Heatmaps"
	if in(Phase, ["SU-Singlet","FakeSU-Singlet","SU-Triplet","FakeSU-Triplet"])
		DirPathOut *= "_Syms=$(Syms...)"
	end
	DirPathOut *= "/xVar=" * xVar * "_yVar=" * yVar * "/"
	mkpath(DirPathOut)

	# Import data
	DF = ImportData(FilePathIn)

	# Unique variables dictionary
	uDict::Dict{String,Vector{Any}} = Dict([
		Var => unique(DF[!, Symbol(Var)]) for Var in AllVars
	])

	xx::Vector{Float64} = uDict[xVar]
	NumX::Int64 = length(xx)
	uDict[xVar] = [NaN]                   # Remove from following cycle
	yy::Vector{Float64} = uDict[yVar]
	NumY::Int64 = length(yy)
	uDict[yVar] = [NaN]                   # Remove from following cycle

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
		for (w,t) in enumerate(uDict["t"]),
			(u,U) in enumerate(uDict["U"]),
			(v,V) in enumerate(uDict["V"]),
			(l,L) in enumerate(uDict["Lx"]),
			(d,δ) in enumerate(uDict["δ"]),
			(b,β) in enumerate(uDict["β"])

			# Initialize local dataframe
			lDict::Dict{String,Any} = Dict([
				"t" => t,
				"U" => U,
				"V" => V,
				"Lx" => L,
				"δ" => δ,
				"β" => β
			])

			# Select entries
			Selections::Vector{Bool} = [true for _ in 1:length(DF.t)]
			for Var in AllVars
				if !in(Var, [xVar, yVar])
					Selections = Selections .* (DF[!, Symbol(Var)] .== lDict[Var])
				end
			end

			# Initialize plot
			H = plot(
				size = (600,400),
				xlabel = L"$%$(VarLabels[xVar])$",
				ylabel = L"$%$(VarLabels[yVar])$",
				legend = :outertopright
			)

			# Write terminal message and file name
			TerminalMsg::String = "\e[2K\e[1GPlotting $(Phase) HF $(HF) data: "
			FilePathOut::String = DirPathOut * "/" * HF
			rawTitle::String = HFLabel * " ("
			for Var in AllVars
				if !in(Var, [xVar, yVar])
					lVar = lDict[Var]
					TerminalMsg *= Var * "=$(lVar), "
					FilePathOut *= "_" * Var * "=$(lVar)"
					RS::String=""
					if !RenormalizeBands && Var=="t"
						RS *= "\\tilde{t}="
					end
					rawTitle *= "\$$(VarLabels[Var])=" * RS * "$(lVar)\$, "
					if Var=="Lx" # I am desperate about correct formatting
						rawTitle = rawTitle[1:end-5] * "\$, "
					elseif Var=="β" && β==Inf
						rawTitle = rawTitle[1:end-6] * "\\infty\$, "
					end
				end
			end
			TerminalMsg = TerminalMsg[1:end-2] * " [x variable: " * xVar *
				", y variable: " * yVar * "]"
			FilePathOut *= "." * Extension
			rawTitle = rawTitle[1:end-2] * ")"
			printstyled("\e[2K\e[1G" * TerminalMsg, color=:yellow)

			# Define x, y variables and h
			vv = DF.v[Selections]
			hh::Matrix{Float64} = zeros(NumY,NumX)
			for j in 1:NumX
				hh[:,j] .= [eval(Meta.parse(
					vv[(j-1) * NumY + i]
				))[HF] for i in 1:NumY]
			end

			# Plot parametrically
			heatmap!(
				xx, yy, abs.(hh),
				color=cs
			)
			title!(L"%$(rawTitle)")

			# Save figure
			Plots.PGFPlotsX.push_preamble!(
				backend_object(H).the_plot, "\\usepackage{bm}"
			)
			savefig(H, FilePathOut)
		end
	end
	printstyled("\e[2K\e[1GDone! Plots saved at $(DirPathOut)\n", color=:green)
end

@doc raw"""
function PlotRMPs(
	Phase::String,
	FilePathIn::String,
	DirPathOut::String;
	xVar::String=\"U\",
	yVar::String=\"V\",
	cs::Symbol=:imola50,
	Extension::String="pdf"
)

Returns: none (plots saved at `DirPathOut`).

`PlotRMPs` (Renormalized Model Parameters) takes as input `Phase` (string
specifying the mean-field phase, the allowed are \"AF\", \"SU-Singlet\",
\"SU-Triplet\"), `FilePathIn` (path to the data files), `DirPathOut` (path to
the output directory). The optional parameters are `xVar` and `yVar` (strings
specifying respectively the x variable and the y variable of the plot, the
allowed are \"t\", \"U\", \"V\", \"Lx\", \"δ\", \"β\"), `cs` (colorscheme
symbol). `Extension` selects the file extension (the allowed are \"pdf\",
\"png\", \"svg\").
"""
function PlotRMPs(
	Phase::String,						# Mean field phase
	FilePathIn::String,					# Data filepath
	DirPathOut::String;					# Output directory path
	xVar::String="U",                   # Specify x variable
	yVar::String="V",                   # Specify y variable
	cs::Symbol=:imola50,                # Custom colorscheme
	Extension::String="pdf"             # File extension
)

	FilePathOut = ""
	AllVars = ["t", "U", "V", "Lx", "δ", "β"]

	# Input safecheck
	if !in(xVar, AllVars)
		@error "Invalid x variable, choose one of $(AllVars)"
		exit()
	end

	if !in(yVar, AllVars)
		@error "Invalid parametric variable, choose one of $(AllVars)"
		exit()
	end

	if xVar==yVar
		@error "You have chosen xVar=yVar!"
		exit()
	end

	# Get LaTeX formatted labels
	VarLabels = GetLabels(Phase)

	# Initialize directory structure
	DirPathOut *= "RMPs"
	if Phase=="SU-Singlet"
		DirPathOut *= "_Syms=$(Syms...)"
	elseif Phase=="SU-Triplet"
		DirPathOut *= "_Syms=$(Syms...)"
	end
	DirPathOut *= "/xVar=" * xVar * "_yVar=" * yVar * "/"
	mkpath(DirPathOut)

	# Import data
	DF = ImportData(FilePathIn)

	# Unique variables dictionary
	uDict::Dict{String,Vector{Any}} = Dict([
		Var => unique(DF[!, Symbol(Var)]) for Var in AllVars
	])

	xx::Vector{Float64} = uDict[xVar]
	NumX::Int64 = length(xx)
	uDict[xVar] = [NaN]                   # Remove from following cycle
	yy::Vector{Float64} = uDict[yVar]
	NumY::Int64 = length(yy)
	uDict[yVar] = [NaN]                   # Remove from following cycle

	# List HF parameters #TODO Generalize
	ListRMPs::Vector{String} = ["reΔ_tilde", "imΔ_tilde", "t_tilde"]

	# Cycle over renormalized model parameters
	for RMP in ListRMPs

		RMPLabel::String = "\$" * VarLabels[RMP] * "\$"

		# Cycle over simulated points
		for (w,t) in enumerate(uDict["t"]),
			(u,U) in enumerate(uDict["U"]),
			(v,V) in enumerate(uDict["V"]),
			(l,L) in enumerate(uDict["Lx"]),
			(d,δ) in enumerate(uDict["δ"]),
			(b,β) in enumerate(uDict["β"])

			# Initialize local dataframe
			lDict::Dict{String,Any} = Dict([
				"t" => t,
				"U" => U,
				"V" => V,
				"Lx" => L,
				"δ" => δ,
				"β" => β
			])

			# Select entries
			Selections::Vector{Bool} = [true for _ in 1:length(DF.t)]
			for Var in AllVars
				if !in(Var, [xVar, yVar])
					Selections = Selections .* (DF[!, Symbol(Var)] .== lDict[Var])
				end
			end

			# Initialize plot
			S = plot(
				size = (600,400),
				 xlabel = L"$%$(VarLabels[xVar])$",
				 ylabel = L"$%$(VarLabels[yVar])$",
				 zlabel = L"$%$(VarLabels[RMP])$",
				legend = :outertopright
			)

			# Write terminal message and file name
			TerminalMsg::String = "\e[2K\e[1GPlotting $(Phase) RMP $(RMP) " *
				"data: "
			FilePathOut::String = DirPathOut * "/" * RMP
			rawTitle::String = RMPLabel * " ("
			for Var in AllVars
				if !in(Var, [xVar, yVar])
					lVar = lDict[Var]
					TerminalMsg *= Var * "=$(lVar), "
					FilePathOut *= "_" * Var * "=$(lVar)"
					RS::String=""
					if !RenormalizeBands && Var=="t"
						RS *= "\\tilde{t}="
					end
					rawTitle *= "\$$(VarLabels[Var])=" * RS * "$(lVar)\$, "
					if Var=="Lx" # I am desperate about correct formatting
						rawTitle = rawTitle[1:end-5] * "\$, "
					elseif Var=="β" && β==Inf
						rawTitle = rawTitle[1:end-6] * "\\infty\$, "
					end
				end
			end
			TerminalMsg = TerminalMsg[1:end-2] * " [x variable: " * xVar *
				", y variable: " * yVar * "]"
			FilePathOut *= "." * Extension
			rawTitle *= " \$k_\\ell=\\pi/3\$)"
			printstyled("\e[2K\e[1G" * TerminalMsg, color=:yellow)

			# Define x, y variables and h
			vv = DF.v[Selections]
			hh::Matrix{Float64} = zeros(NumY,NumX)

			# Plot individually in order to adjust visual angles
			if RMP=="reΔ_tilde"

				for j in 1:NumX
					hh[:,j] .= [eval(Meta.parse(
						vv[(j-1) * NumY + i]
					))["m"] for i in 1:NumY]
				end

				# Plot parametrically
				zz = hh .* (xx' .+ 8*yy)
				surface!(
					xx, yy, zz,
					color=cs,
					label=L"$%$(VarLabels[RMP])$",
					camera=(30,25),
#                    zlim=(0.65,1),
#                    clim=(0.65,1)
				)
				title!(L"%$(rawTitle)")

			elseif RMP=="imΔ_tilde"

				for j in 1:NumX
					hh[:,j] .= [eval(Meta.parse(
						vv[(j-1) * NumY + i]
					))["wp"] for i in 1:NumY]
				end

				# Plot parametrically
				zz = hh .* (2 .* yy * ones(1,length(xx)))
				surface!(
					xx, yy, zz,
					color=cs,
					label=L"$%$(VarLabels[RMP])$",
					camera=(30,25),
#                    zlim=(0.65,1),
#                    clim=(0.65,1)
				)
				title!(L"%$(rawTitle)")

			elseif RMP=="t_tilde"

				for j in 1:NumX
					hh[:,j] .= [eval(Meta.parse(
						vv[(j-1) * NumY + i]
					))["w0"] for i in 1:NumY]
				end

				# Plot parametrically
				zV::Matrix{Float64} = zeros(size(hh))
				if xVar=="V"
					zV = ones(length(yy),1) * xx'
				elseif yVar=="V"
					zV = yy * ones(1,length(xx))
				end
				zz = t*ones(size(hh)) .- hh .* zV
				surface!(
					xx, yy, zz,
					color=cs,
					label=L"$%$(VarLabels[RMP])$",
					camera=(80,25),
					zlim=(0.75,1),
					clim=(0.75,1)
				)
				title!(L"%$(rawTitle)")
			end

			# Save figure
			Plots.PGFPlotsX.push_preamble!(
				backend_object(S).the_plot, "\\usepackage{bm}"
			)
			savefig(S, FilePathOut)
		end
	end
	printstyled("\e[2K\e[1GDone! Plots saved at $(DirPathOut)\n", color=:green)
end