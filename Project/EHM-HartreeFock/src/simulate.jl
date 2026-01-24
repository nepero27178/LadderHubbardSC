#!/usr/bin/julia
using DelimitedFiles
using Dates
using DataFrames

# Arguments handler
if length(ARGS) != 1
	println("How to use this program?
Type the following: \$ julia simulate.jl --mode
Where:
· mode = hs / rs")
	exit()
else
	UserInput = ARGS
	Mode = UserInput[1][3:end]
end

# Includer
PROJECT_SRC_DIR = @__DIR__
if in(Mode, ["rs", "hs"])
	include(PROJECT_SRC_DIR * "/setup/" * Mode * "-setup.jl")
else
	@error "Invalid argument. Use: mode = hs / rs"
end
include(PROJECT_SRC_DIR * "/modules/methods-simulating.jl")
include(PROJECT_SRC_DIR * "/modules/methods-physics.jl")
include(PROJECT_SRC_DIR * "/modules/methods-optimizations.jl")
include(PROJECT_SRC_DIR * "/modules/methods-IO.jl")

# Routines
@doc raw"""
function RunHFScan(
	UU::Vector{Float64},
	VV::Vector{Float64},
	LL::Vector{Int64},
	δδ::Vector{Float64},
	ββ::Vector{Float64},
	p::Int64,
	Δm::Dict{String,Float64},
	Δn::Float64,
	g::Float64;
	FilePathOut::String="",
	RenormalizeBands::Bool=true
)

Returns: none if `FilePathOut` is specified.

`RunHFScan` takes as input `UU` (vector of local repulsions), `VV` (vector of
non local attractions), `LL` (vector of the square lattice dimensions), `δδ` 
(vector of dopings with respect to the half filled lattice), `ββ` (vector of 
inverse temperatures), `p` (maximum number of HF iterations), `Δm` (tolerance on
each order parameter) and `Δn` (tolerance on density in chemical potential
estimation), `g` (mixing parameter). It performs an iterative HF analysis over a
sequence of 2D square lattices for all the possible combinations of the 
specified parameters. Check the source files in the `/modules` folder for more
informations on the algorithm. The boolean option `RenormalizeBands' allows
for choosing to renormalize or not the hopping parameter.
"""
function RunHFScan(
	Phase::String,						# Mean field phase
	tt::Vector{Float64},				# Hopping amplitude
	UU::Vector{Float64},				# Local repulsion
	VV::Vector{Float64},				# Non-local attraction
	LL::Vector{Int64},					# Lattice size
	δδ::Vector{Float64},				# Doping
	ββ::Vector{Float64},				# Inverse temperature
	p::Int64,							# Number of iterations
	Δv::Dict{String,Float64},			# Tolerance on magnetization
	Δn::Float64,						# Tolerance on density
	g0::Float64;						# Mixing parameter
	Syms::Vector{String}=["s"],			# Gap function symmetries
	FilePathOut::String="",				# Output file
	InitializeFile::Bool=true,			# Initialize file at FilePathOut
	RenormalizeBands::Bool=true,		# Conditional renormalization of t
	OptimizeBZ::Bool=true,				# Conditional optimization of BZ
	Optimizeg::Bool=true				# Conditional optimization of g
)

	# Warn user of memory-heavy simulations detection
	Iterations = length(UU) * length(VV) * length(LL) * length(ββ) * length(δδ)
	if FilePathOut == "" && Iterations > 200
		@warn "No output file specified and more than 200 simulations " *
			"request detected. Simulations results are going to be stored " *
			"in your memory. Consider specifying a `FilePathOut` and " *
			"unloading your memory."
	end

	# Get Hartree Fock Parameters labels
	HFPs = GetHFPs(Phase;Syms)

	# File coditional initialization (otherwise, just append)
	if FilePathOut != "" && InitializeFile
		Header = "t;U;V;Lx;β;δ;v;Q;ΔT;I;μ;g0;g;f;f0\n"
		write(FilePathOut, Header)
	end

	# Initializers
	v0::Dict{String,Float64} = Dict([
		key => 0.1 for key in HFPs
	])

	# HF iterations
	i = 1
	for t in tt,
		Lx in LL,
		δ in δδ,
		β in ββ

		Uc::Float64 = 0.0
		if Optimizeg && occursin("SU", Phase)
			Uc = GetUc(t,[Lx,Lx],δ,β)
		end

		g = g0
		for U in UU # TODO Evaluate possibility of computing dynamically g with the bare bands

			if U>(2/g0-1)*Uc && Optimizeg && occursin("SU", Phase)
				Og = GetOptimalg(U,Uc)
				Og<g0 ? g=Og : 0
			end

			for V in VV
			#NOTE Due to order inversion here (Vδ) => (δV), heatmap now gets: xVar=δ, yVar=V

				Parameters::Dict{String,Float64} = Dict([
					"t" => t,
					"U" => U,
					"V" => V
				])

				L = [Lx, Lx]
				printstyled(
					"\e[2K\e[1GRun ($i/$Iterations): " *
					"$Phase HF at t=$t, U=$U, V=$V, L=$Lx, β=$β, δ=$δ",
					color=:yellow
				)
				ResultsVector::Matrix{Any} = [0 0]  # Dummy placeholder
				ResultsVector = hcat(ResultsVector, [t U V 0 β δ]) # Lx later

				# Run routine, all positional arguments here must be false
				Run, Performance = RunHFAlgorithm(
					Phase,Parameters,L,0.5+δ,β,
					p,Δv,Δn,g;
					v0i=v0,
					Syms,
					RenormalizeBands,
					OptimizeBZ
				)

				v::Dict{String,Float64} = Dict([
					key => Run["HFPs"][key] for key in HFPs
				])

				fMFT = Run["FreeEnergy"]
				Qs::Dict{String,Float64} = Dict([
					key => Performance["Quality"][key] for key in HFPs
				])
				ResultsVector = hcat(
					ResultsVector[:,3:end],
					[v Qs Performance["Runtime"] Performance["Steps"] Run["ChemicalPotential"] g0 g fMFT]
				)
				ResultsVector[4] = Lx # Add here to get correct Int64 formatting

				i += 1

				# Append to initialized or existing file
				if FilePathOut != ""
					open(FilePathOut, "a") do io
						writedlm(io, ResultsVector, ';')
					end
				end

				# Use current step as initializer coherently (no NaN and no 0)
				if any(isnan.(values(v)) .|| values(v).==0.0)
					0 # println(v0)
				else
					v0 = copy(v)
				end
			end
		end
	end

	printstyled(
		"\e[2K\e[1GDone! Data saved at " * FilePathOut * "\n", color=:green
	)

end

# Main run
function main()

	# Create output directory
	# For simulations: Setup > Phase > Syms (to make comparable data in the same folder)
	# For plots: Phase > Setup > Syms (to make same-phase plots in the same folder)
	DirPathOut = PROJECT_SRC_DIR * "/../simulations/Mode=$(Mode)/Setup=$(Setup)/Phase=$(Phase)"
	FilePathOut = DirPathOut * "/Syms=$(Syms...).csv"
	mkpath(dirname(FilePathOut))

	# Filter out non half-filled simulations from AF phase
	occursin("AF", Phase) ? filter!(==(0),δδ) : 0

	TotalRunTime = @elapsed begin
		RunHFScan(
			Phase,
			tt,UU,VV,
			LL,δδ,ββ,
			p,Δv,Δn,g;
			Syms,
			FilePathOut,
			RenormalizeBands,
			OptimizeBZ=false # For now
		)
	end

	LogPathOut = DirPathOut * "/Syms=$(Syms...)_log.txt"
	Header = "tt;UU;VV;LL;ββ;δδ;p;Δv;Δn;g;TotalRunTime;Machine\n"
	write(LogPathOut, Header)
	Log = [[tt] [UU] [VV] [LL] [ββ] [δδ] p Δv Δn TotalRunTime gethostname()]
	open(LogPathOut, "a") do io
		writedlm(io, Log, ';')
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
