# Includer
PROJECT_ROOT = @__DIR__
PROJECT_ROOT *= "/.."   # Up to the effective root
include(PROJECT_ROOT * "/setup/graphic-setup.jl")
include(PROJECT_ROOT * "/modules/methods-simulating.jl")

@doc raw"""
[TODO]
"""
function PlotDensity(
		Phase::String,
		U::Float64,
		Δ::Float64,
		Lx::Int64,
		ββ::Vector{Float64},
		FilePathOut::String;
		cs::Symbol=:winter,
	)

	P = plot(
		size = (600,400),
		xlabel = L"$\mu/\Delta$",
		ylabel = L"$N(\mu)/2L_xL_y$",
		xticks = [-4,-2,-1,0,1,2,4],
		xminorticks = false,
		ylim = (-0.05,1.05),
		title = L"%$(Phase) density ($U=%$(U)$, $L=%$(Lx)$)",
		legend = :topleft,
		legendfonthalign = :left
	)

	Parameters::Dict{String,Float64} = Dict([
		"t" => 1.0,
		"U" => U
	])

	D = 2*Lx^2

	# Reciprocal space discretization (normalized to 1)
	Kx::Vector{Float64} = [kx for kx in -1:2/Lx:1]
	popfirst!(Kx)
	Ky::Vector{Float64} = [ky for ky in -1:2/Lx:1]
	popfirst!(Ky)
	K::Matrix{Vector{Float64}} = [ [kx,ky] for kx in Kx, ky in Ky ]

	v::Dict{String,Float64} = Dict("m" => Δ/U)
	vline!(
		[-Δ,Δ],
		ls = :dash,
		color = :gray,
		alpha = 0.5,
		label = "",
	)
	annotate!([
		(-Δ-0.25,0.75,text(L"\mu=-\Delta",9,color=:gray,alpha=0.5,rotation=90)),
		(Δ+0.25,0.25,text(L"\mu=\Delta",9,color=:gray,alpha=0.5,rotation=90))
	])

	xx = [x for x in -5:0.01:5]
	q = floor(Int64, length(colorschemes[cs]) / length(ββ) )
	for (j,β) in enumerate(ββ)
		n(μ::Float64) = sum( GetKPopulation(Phase,Parameters,K,v,μ,β) )/D
		yy = n.(xx)
		if β==Inf
			βlabel = "\\beta=\\infty"
		else
			βlabel = "\\beta=$(β)"
		end
		plot!(
			xx,yy,
			label=L"%$(βlabel)",
			color=colorschemes[cs][q*j]
		)
	end

	# Save figure
	savefig(P, FilePathOut)
end

if abspath(PROGRAM_FILE) == @__FILE__
	Phase = "AF"
	U = 10.0
	Δ = 1.0
	Lx = 256
	ββ = [Inf, 100.0, 10.0, 1.0, 0.1]
	FilePathOut = "n_$(Phase)_U=$(U)_Δ=$(Δ).pdf"
	PlotDensity(Phase,U,Δ,Lx,ββ,FilePathOut)
end