using DataFrames
using DelimitedFiles

PROJECT_LABS_DIR = @__DIR__
include(PROJECT_LABS_DIR * "/../../modules/methods-simulating.jl")

@doc raw"""
function FindFS(
	Phase::String,
	Parameters::Dict{String,Float64},
	K::Matrix{Vector{Float64}},
	v::Dict{String,Float64},
	μ::Float64;
	RenormalizeBands::Bool=true
)::Tuple{Vector{Float64}, Vector{Float64}}

Returns: list of `x`, `y` coordinates of the Fermi Surface.

[...]
"""
function FindFS(
	Phase::String,						# Mean field phase
	Parameters::Dict{String,Float64},	# Model parameters t,U,V
	K::Matrix{Vector{Float64}},			# BZ grid
	v::Dict{String,Float64},				# HF parameters
	μ::Float64;							# Chemical potential
	RenormalizeBands::Bool=true		# Conditional renormalization of t
)::Tuple{Vector{Float64}, Vector{Float64}}

	# Initialize collectors
	xx = Float64[]
	yy = Float64[]

	t = Parameters["t"]
	if RenormalizeBands
		# Conditional renormalization of bands
		t -= v["w0"] * Parameters["V"]
	end

	reΔk::Float64 = v["m"] * (Parameters["U"] + 8*Parameters["V"])
	Q = K .* pi
	for i in 1:(size(Q,1)-1), j in 1:(size(Q,2)-1)

		# Check all four edges of the cell
		for ((x1, y1), (x2, y2)) in [
			(Q[j,i], Q[j,i+1]),
			(Q[j,i+1], Q[j+1,i+1]),
			(Q[j+1,i+1], Q[j+1,i]),
			(Q[j+1,i], Q[j,i]),
		]

			if in(Phase, ["AF", "FakeAF"])

				δE::Vector{Float64} = zeros(2)
				for (i,k) in enumerate([[x1, y1], [x2, y2]])
					# Renormalized bands
					εk::Float64 = GetHoppingEnergy(t,k)

					# Renormalized gap
					imΔk::Float64 = 2 * v["wp"] * Parameters["V"] * StructureFactor("S",k)

					# Renormalized gapped bands
					Ek::Float64 = sqrt( εk^2 + reΔk^2 + imΔk^2 )

					δE[i] = Ek - μ
				end

				# If sign alternates: estimate crossing
				if prod(δE) < 0

					# Parametric position on the edge
					k(g) = [x1 + g*(x2-x1), y1 + g*(y2-y1)]

					# Functional to be minimized
					F(g) = Float64(sqrt(
							( GetHoppingEnergy(t,k(g)) )^2 +
							( reΔk )^2 +
							( 2 * v["wp"] * Parameters["V"] * StructureFactor("S",k(g)) )^2
						)) - μ

					# Bisect until convergence
					g = find_zero(F, (0, 1), Bisection())
					push!(xx, x1 + g*(x2-x1))
					push!(yy, y1 + g*(y2-y1))

				end
			end
		end
	end

	return xx, yy
end

@doc raw"""
function GetPointFS(
	Phase::String,
	t::Float64,
	U::Float64,
	V::Float64,
	Lx::Int64,
	δ::Float64,
	β::Float64,
	FilePathIn::String,
	DirPathOut::String
)::DataFrame

Returns: full information on point and relative Fermi Surface.

[...]
"""
function GetPointFS(
	Phase::String,						# Mean field phase
	t::Float64,							# Hopping amplitude
	U::Float64,							# Local repulsion
	V::Float64,							# Non-local attraction
	Lx::Int64,							# Lattice size
	δ::Float64,							# Doping
	β::Float64,							# Inverse temperature
	FilePathIn::String,					# Data filepath
	DirPathOut::String;					# Output directory path
)::DataFrame

	# Initialize directory structure
	DirPathOut *= "/Lx=$(Lx)_δ=$(δ)_β=$(β)/"
	mkpath(DirPathOut)
	FilePathOut = DirPathOut * "t=$(t)_U=$(U)_V=$(V).txt"

	# Get Hartree Fock Parameters labels
	HFPs = GetHFPs(Phase)

	# Initialize model parameters
	Parameters::Dict{String,Float64} = Dict([
		"t" => t,
		"U" => U,
		"V" => V
	])

	# Read data and create DataFrame
	DataIn = open(FilePathIn) do io
		readdlm(FilePathIn, ',', comments=false, '\n')
	end
	DF = DataFrame(DataIn[2:end,:], DataIn[1,:])
	DF.Lx = Int64.(DF.Lx)
	Point = DF[
			(DF.t.==t) .&
			(DF.U.==U) .&
			(DF.V.==V) .&
			(DF.Lx.==Lx) .&
			(DF.δ.==δ) .&
			(DF.β.==β),
		:]

	# Reciprocal space discretization (normalized to 1)
	Kx::Vector{Float64} = [kx for kx in -1:2/Lx:1]
	popfirst!(Kx)
	Ky::Vector{Float64} = [ky for ky in -1:2/Lx:1]
	popfirst!(Ky)
	K::Matrix{Vector{Float64}} = [ [kx,ky] for kx in Kx, ky in Ky ]

	v = eval.(Meta.parse.(Point.v))[1]
	μ = Point.μ[1]

	RenormalizeBands::Bool = true
	if Phase=="FakeAF"
		RenormalizeBands = false
	end
	xx, yy = FindFS(Phase,Parameters,K,v,μ;RenormalizeBands)
	Point[!, :Phase] .= [Phase]
	Point[!, :Curve] .= [Dict(["kx" => xx, "ky" => yy])]

	writedlm(FilePathOut, [xx yy])
	return Point

end

function PlotPointFS!(
	Phase::String,						# Mean field phase
	t::Float64,							# Hopping amplitude
	U::Float64,							# Local repulsion
	V::Float64,							# Non-local attraction
	DirPathOut::String;					# Output directory path
	FilePathIn::String=""
)
	xx = Float64[]
	yy = Float64[]
	if !any(readdir(DirPathOut) .== "$(Phase)_t=$(t)_U=$(U)_V=$(V).txt")
		PointFS = GetPointFS(
			"AF",
			t,U,V,
			Lx,δ,β,
			FilePathIn,
			DirPathOut
		)
		xx = PointFS.Curve[1]["kx"]./pi
		yy = PointFS.Curve[1]["ky"]./pi

	elseif any(readdir(DirPathOut) .== "$(Phase)_t=$(t)_U=$(U)_V=$(V).txt")
		FilePathIn = DirPathOut * "$(Phase)_t=$(t)_U=$(U)_V=$(V).txt"
		Data = open(FilePathIn) do io
			readdlm(FilePathIn, '\t', comments=false, '\n')
		end
		xx = Data[:,1]
		yy = Data[:,2]
	end

	scatter!(
		xx,yy,
		markerstrokealpha = 0,
		markersize = 1.5,
		label=Phase
	)

	return xx, yy

end

if abspath(PROGRAM_FILE) == @__FILE__
	Setup = "B[256]"
	DirPathOut = PROJECT_ROOT * "/../processing/fermi-surface/"

	Lx = 256
	δ = 0.3
	β = 100.0
	DirPathOut *= "/Lx=$(Lx)_δ=$(δ)_β=$(β)/"

	t = 1.0
	U = 20.0
	V = 3.0

	S = scatter(
		markerstrokealpha = 0,
		markersize = 1.5,
		size = (500,500),
		aspect_ratio = :equal,
		legendfonthalign = :left,
	)

	for Phase in ["AF", "FakeAF"]
		FilePathIn = PROJECT_ROOT * "/simulations/scan/Setup=$(Setup)/$(Phase).txt"
		PlotPointFS!(Phase,t,U,V,DirPathOut;FilePathIn)
	end

	title!(L"U=%$(U), V=%$(V), \delta=%$(δ)")
	xlabel!(L"k_x/\pi")
	ylabel!(L"k_y/\pi")
	savefig(S, DirPathOut * "t=$(t)_U=$(U)_V=$(V).pdf")
end
