# --- Fermi Surface finder ---

@doc raw"""
function FindFS(
	Phase::String,
	Parameters::Dict{String,Float64},
	K::Matrix{Vector{Float64}},
	v::Dict{String,Float64},
	μ::Float64;
	RenormalizeHopping::Bool=true
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
	RenormalizeHopping::Bool=true		# Conditional renormalization of t
)::Tuple{Vector{Float64}, Vector{Float64}}

	# Initialize collectors
	xx = Float64[]
	yy = Float64[]

	t = Parameters["t"]
	if RenormalizeHopping
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

	RenormalizeHopping::Bool = true
	if Phase=="FakeAF"
		RenormalizeHopping = false
	end
	xx, yy = FindFS(Phase,Parameters,K,v,μ;RenormalizeHopping)
	Point[!, :Phase] .= [Phase]
	Point[!, :Curve] .= [Dict(["kx" => xx, "ky" => yy])]

	writedlm(FilePathOut, [xx yy])
	return Point

end

# --- Post-computation of μ ---

@doc raw"""
function Addμ(
	Phase::String,
	FilePathIn::String,
	FilePathOut::String
)::DataFrame

Returns: new dataframes completed with chemical potential.

[...]
"""
function Addμ(
	Phase::String,
	FilePathIn::String,
	FilePathOut::String
)::DataFrame

	@warn "This is a temporary function and is going to be suppressed."

	# Read data and create DataFrame
	DataIn = open(FilePathIn) do io
		readdlm(FilePathIn, ';', comments=false, '\n')
	end
	DF = DataFrame(DataIn[2:end,:], DataIn[1,:])
	DF.Lx = Int64.(DF.Lx)

	μμ::Vector{Float64} = zeros(length(DF.t))
	RenormalizeHopping::Bool = true
	if Phase=="FakeAF"
		RenormalizeHopping = false
	end

	for i in 1:length(μμ)

		# Model parameters
		Parameters::Dict{String,Float64} = Dict([
			"t" => DF.t[i],
			"U" => DF.U[i],
			"V" => DF.V[i]
		])

		Lx = DF.Lx[i]
		# Reciprocal space discretization (normalized to 1)
		Kx::Vector{Float64} = [kx for kx in -1:2/Lx:1]
		popfirst!(Kx)
		Ky::Vector{Float64} = [ky for ky in -1:2/Lx:1]
		popfirst!(Ky)
		K::Matrix{Vector{Float64}} = [ [kx,ky] for kx in Kx, ky in Ky ]

		β = DF.β[i]
		v = eval(Meta.parse(DF.v[i]))
		n = 0.5 + DF.δ[i]

		μμ[i] = FindRootμ(Phase,Parameters,K,v,n,β;Δn=1e-2,RenormalizeHopping)
	end

	DF[!,:μ] = μμ
	writedlm(FilePathOut, Iterators.flatten(([names(DF)],eachrow(DF))), ';')

	return DF
end