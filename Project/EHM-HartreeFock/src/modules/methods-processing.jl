PROJECT_METHODS_DIR = @__DIR__
include(PROJECT_METHODS_DIR * "/methods-simulating.jl")

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

@doc raw"""
function Addμ(
	Phase::String,
	FilePathIn::String,
	FilePathOut::String
)::DataFrame

Returns: new dataframes completed with chemical potential.
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
	RenormalizeBands::Bool = true
	if Phase=="FakeAF"
		RenormalizeBands = false
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

		μμ[i] = FindRootμ(Phase,Parameters,K,v,n,β;Δn=1e-2,RenormalizeBands)
	end

	DF[!,:μ] = μμ
	writedlm(FilePathOut, Iterators.flatten(([names(DF)],eachrow(DF))), ';')

	return DF
end

function FilteredRun(
	FilePathIn::String
)

	# Unpack path, load data and log
	Setup, Phase, Syms = UnpackFilePath(FilePathIn)
	DF = ImportData(FilePathIn)
	LogPathIn = replace(FilePathIn, ".csv" => ".log")
	LogDF = ImportData(LogPathIn)

	# Get Hartree Fock Parameters labels
	HFPs = GetHFPs(Phase;Syms)

	# Set up bands renormalization
	RenormalizeBands::Bool = true
	occursin("Fake",Phase) ? RenormalizeBands = false : 0

	# Prepare for filtered run
	Δv::Dict{String,Float64} = eval(Meta.parse(LogDF.Δv[1]))
	Δn::Float64 = LogDF.Δn[1]
	p::Int64 = 3*LogDF.p[1] # Four times as many steps

	# Build output
	FilePathOut = replace(FilePathIn, "simulations" => "processing")
	LogPathOut = replace(LogPathIn, "simulations" => "processing")
	mkpath(dirname(FilePathOut))

	# Assert if working on an already filtered layer
	if occursin("Layer=",FilePathIn)
		eval(Meta.parse(
			split(FilePathIn,'.')[2]
		))
		FilePathOut = replace(FilePathOut,"Layer=$(Layer)" => "Layer=$(Layer+1)")
		LogPathOut = replace(LogPathOut,"Layer=$(Layer)" => "Layer=$(Layer+1)")
	elseif !occursin("Layer=",FilePathIn)
		FilePathOut = replace(FilePathOut,".csv" => ".Layer=1.csv")
		LogPathOut = replace(LogPathOut,".log" => ".Layer=1.log")
	end

	# Write output header
	Header = "t;U;V;Lx;β;δ;v;Q;ΔT;I;μ;g0;g;fMFT\n"
	write(FilePathOut, Header)

	# Write log
	Header = "p;Δv;Δn;Machine\n"
	write(LogPathOut, Header)
	Log = [p Δv Δn gethostname()]
	open(LogPathOut, "a") do io
		writedlm(io, Log, ';')
	end

	j::Int64 = 1 # NaN counter
	J::Int64 = 0 # Remaining NaN counter
	Iterations::Int64 = sum(isnan.(DF.fMFT))
	i0 = findfirst(!isnan,DF.fMFT) # Hopefully i0=1
	v0i::Dict{String,Float64} = eval(Meta.parse(DF[i0,:].v))
	for (r,row) in enumerate(eachrow(DF))

		if isnan(row.fMFT) # Check energy to assess convergence
			printstyled(
				"\e[2K\e[1GFilter run ($j/$Iterations): " *
				"$Phase HF at t=$(row.t), U=$(row.U), V=$(row.V), L=$(row.Lx), β=$(row.β), δ=$(row.δ)",
				color=:yellow
			)
			ResultsVector::Matrix{Any} = [0 0]  # Dummy placeholder
			ResultsVector = hcat(ResultsVector, [row.t row.U row.V 0 row.β row.δ]) # Lx later

			Parameters::Dict{String,Float64} = Dict("t" => row.t, "U" => row.U, "V" => row.V)
			for (key,val) in v0i
				v0i[key] = val + 0.1*sign(val) # Avoid 0.0 prolems
			end

			g::Float64 = row.g * 0.75
			Run, Performance = RunHFAlgorithm(
				Phase,Parameters,[row.Lx, row.Lx],0.5+row.δ,Float64(row.β),
				p,Δv,Δn,g;
				v0i,
				Syms,
				RenormalizeBands,
				OptimizeBZ=false # More precise
			)

			v::Dict{String,Float64} = Dict([
				key => Run["HFPs"][key] for key in HFPs
			])
			fMFT::Float64 = Run["FreeEnergy"]
			Qs::Dict{String,Float64} = Dict([
				key => Performance["Quality"][key] for key in HFPs
			])
			ResultsVector = hcat(
				ResultsVector[:,3:end],
				[v Qs Performance["Runtime"] Performance["Steps"] Run["ChemicalPotential"] row.g0 g fMFT]
			)
			ResultsVector[4] = row.Lx # Add here to get correct Int64 formatting

			open(FilePathOut, "a") do io
				writedlm(io, ResultsVector, ';')
			end

			isnan(fMFT) ? J += 1 : v0i = copy(v) # If converged, advance
			j += 1

		elseif !isnan(row.fMFT) # Check energy for convergence

			ResultsVector = [0 0]  # Dummy placeholder
			ResultsVector = hcat(ResultsVector, [row.t row.U row.V 0 row.β row.δ]) # Lx later
			v = eval(Meta.parse(DF[r,:].v))
			v0i = copy(v)
			Qs = eval(Meta.parse(DF[r,:].v))
			ResultsVector = hcat(
				ResultsVector[:,3:end],
				[v Qs row.ΔT row.I row.μ row.g0 row.g row.fMFT]
			)
			ResultsVector[4] = row.Lx # Add here to get correct Int64 formatting

			open(FilePathOut, "a") do io
				writedlm(io, ResultsVector, ';')
			end
		end
	end

	printstyled(
		"\e[2K\e[1GDone! Data saved at " * FilePathOut * "\n", color=:green
	)
	if J>0
		@warn "Remaining NaNs:" J
	end

end