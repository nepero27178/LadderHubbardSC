using DataFrames
using DelimitedFiles

PROJECT_METHODS_DIR = @__DIR__
include(PROJECT_METHODS_DIR * "/structs.jl")
include(PROJECT_METHODS_DIR * "/methods-physics.jl")

function ImportData(
    FilePathIn::String
)::DataFrame

	# Read data and create DataFrame
	DataIn = open(FilePathIn) do io
		readdlm(FilePathIn, ';', comments=false, '\n')
	end
	DF = DataFrame(DataIn[2:end,:], DataIn[1,:])
    DF = identity.(DF) # Format columns type
    # DF.Lx = Int64.(DF.Lx)

	return DF
end

	
function UnpackFilePath(
	FilePathIn::String
)::Tuple{String,String,Vector{String}}

	SubStrs = split(FilePathIn,'/')
	SetupStr = SubStrs[end-2]
	PhaseStr = SubStrs[end-1]
	SymsStr = SubStrs[end]
	
	Setup = split(SetupStr,'=')[2]
	Phase = split(PhaseStr,'=')[2]
	Syms = split(split(SymsStr,'=')[2],'.')[1]

	return String(Setup), String(Phase), [string(s) for s in Syms]
end
	
function ReshapeData(
	DF::DataFrame;
	xVar::String="U",
	yVar::String="V",
	zVar::String="fMFT"
)::Tuple{Any, Any, Any}
	
	xx = unique(DF[!,xVar])
	yy = unique(DF[!,yVar])
	zz = reshape(DF[!,zVar],length(yy),length(xx))
	
	return xx, yy, zz
end

function EnlargeDF!(
	Sim::Simulation
)::DataFrame
	
	DF = Sim.DF

	# List HF parameters
	ListHFPs::Vector{String} = [key for key in keys(
		eval(Meta.parse(
			DF.v[1]
		))
	)]
	vv = eval.(Meta.parse.(DF.v))
	QQ = eval.(Meta.parse.(DF.Q))
	for HF in ListHFPs
		DF[!,HF] = get.(vv,HF,"N/A")
		DF[!,"Q" * HF] = get.(QQ,HF,"N/A")
	end
	
	for (RMP,F) in GetRMPs(Sim.Phase;Sim.Syms)
		DF[!,RMP] = F(DF)
	end

	return select!(DF, Not(:v, :Q))
	
end

@doc raw"""
function GetHFPs(
	Phase::String;
	Syms::Vector{String}=["s"]
)::Vector{String}

Returns: Hartree Fock Parameters labels for the given Phase.
"""
function GetHFPs(
	Phase::String;						# Mean field phase
	Syms::Vector{String}=["s"]			# Gap function symmetries
)::Vector{String}

	AF = false
	Singlet = false
	Triplet = false
	SymErr = "Invalid symmetries. $(Syms) is incoherent with $(Phase)."
	if in(Phase, ["AF", "FakeAF"])
		AF = true
	elseif in(Phase, ["SU-Singlet", "FakeSU-Singlet"])
		issubset(Syms, ["s", "S", "d"]) ? Singlet = true : throw(SymErr)
	elseif in(Phase, ["SU-Triplet", "FakeSU-Triplet"])
		issubset(Syms, ["px", "py", "p+", "p-"]) ? Triplet = true : throw(SymErr)
	end

	KeysList::Vector{String} = []
	AF ? KeysList = ["m", "w0", "wp"] : 0
	Singlet ? KeysList = vcat(["Δ$(Sym)" for Sym in Syms], "gS", "gd") : 0
	Triplet ? KeysList = vcat(["Δ$(Sym)" for Sym in Syms], "gS", "gd") : 0

	# Pure symmetry drop TODO Extension to Triplet
	if Singlet && (sort(Syms) == ["S", "s"] || Syms == ["d"])
		pop!(KeysList) # Pop last: gd
	end

	return KeysList

end

@doc raw"""
function GetRMPs(
	Phase::String,
	Syms::Vector{String}=["s"]
)::Vector{String}

Returns: Renormalized Model Parameters labels for the given Phase.
"""
function GetRMPs(
	Phase::String;						# Mean field phase
	Syms::Vector{String}=["s"]			# Gap function symmetries
)::Dict{String,Any}

	AF = false
	Singlet = false
	Triplet = false
	SymErr = "Invalid symmetries. $(Syms) is incoherent with $(Phase)."
	if in(Phase, ["AF", "FakeAF"])
		AF = true
	elseif Phase=="SU-Singlet"
		issubset(Syms, ["s", "S", "d"]) ? Singlet = true : throw(SymErr)
	elseif Phase=="SU-Triplet"
		issubset(Syms, ["px", "py", "p+", "p-"]) ? Triplet = true : throw(SymErr)
	end

	AF ? KeysList = ["reΔ_tilde", "imΔ_tilde", "t_tilde"] : 0
	Singlet || Triplet ? KeysList = ["t_tilde"] : 0 #TODO Add bands renormalization
	if AF
		# Hopping renormalization
		RtAF(df) = df.t .- df.w0 .* df.V
		RReΔ(df) = df.m .* (df.U .+ 8*df.V)
		RImΔ(df) = 2 * df.wp .* df.V
		RMPs = Dict(
			"Rt" => df -> RtAF(df),
			"RReΔ" => df -> RReΔ(df),
			"RImΔ" => df -> RImΔ(df),
		)
		return RMPs
	elseif Singlet || Triplet
		# Hopping renormalization
		RtSU(df) = df.t .- df.gS/2 .* df.V
		RMPs = Dict(
			"Rt" => df -> RtSU(df)
		)
		return RMPs
	end

end

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
		"g0" => "g_0",
		"g" => "g",
		"fMFT" => "f_\\mathrm{MFT}"
	])

	if in(Phase, ["AF", "FakeAF"])
		PhaseLabels::Dict{String,String} = Dict([
			# HFPs
			"m" => "m",
			"w0" => "w^{(\\mathbf{0})}",
			"wp" => "w^{(\\pi)}",
			# Convergence
			"Qm" => "Q(m)",
			"Qw0" => "Q(w^{(\\mathbf{0})})",
			"Qwp" => "Q(w^{(\\pi)})",
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
			"Rt" => "\\tilde{t}"
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
			"Rt" => "\\tilde{t}"
		])
	end

	return merge(VarLabels, PhaseLabels)

end