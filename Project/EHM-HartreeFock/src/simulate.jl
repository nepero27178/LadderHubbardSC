#!/usr/bin/julia
using DelimitedFiles
using Dates
using DataFrames

# Arguments handler
if length(ARGS) != 1
	println("How to use this program?
Type the following: \$ julia ./simulate.jl --mode
Where:
· --mode = --scan / --heatmap / --record-g")
	exit()
else
	UserInput = ARGS
	Mode = UserInput[1][3:end]
end

# Includer
PROJECT_ROOT = @__DIR__
PROJECT_ROOT *= "/.."   # Up to the effective root
if in(Mode, ["scan", "heatmap", "record-g"])
	include(PROJECT_ROOT * "/src/setup/" * Mode * "-simulations-setup.jl")
else
	@error "Invalid argument. Use: --mode = --scan / --heatmap / --record-g"
end
include(PROJECT_ROOT * "/src/modules/methods-simulating.jl")
include(PROJECT_ROOT * "/src/modules/methods-processing.jl")

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
		Header = "t;U;V;Lx;β;δ;v;Q;ΔT;I;μ;g0;g;fMFT\n"
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

		if Optimizeg && occursin("SU", Phase)
			Uc = GetUc(t,[Lx,Lx],δ,β)
		end

		g = g0
		for U in UU

			if U>(2/g0-1)*Uc && Optimizeg && occursin("SU", Phase)
				Og = GetOptimalg(U,Uc)
				Og<g0 ? g=Og : 0
			end

			for V in VV

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
				ResultsVector = hcat(ResultsVector, [t U V Lx β δ])

				# Run routine, all positional arguments here must be false
				Run, Performance = RunHFAlgorithm(
					Phase,Parameters,L,0.5+δ,β,
					p,Δv,Δn,g;
					# v0i=v0, # TODO Processing: find and replace suspect points
					Syms,
					RenormalizeBands
				)

				v::Dict{String,Float64} = Dict([
					key => Run["HFPs"][key] for key in HFPs
				])

				fMFT = Run["FreeEnergy"]
				Qs::Dict{String,Float64} = Dict([
					key => Performance["Quality"][key] for key in HFPs
				])
				ResultsVector = hcat(ResultsVector[:,3:end], [v Qs Performance["Runtime"] Performance["Steps"] Run["ChemicalPotential"] g0 g fMFT])

				i += 1

				# Append to initialized or existing file
				if FilePathOut != ""
					open(FilePathOut, "a") do io
						writedlm(io, ResultsVector, ';')
					end
				end

				v0 = copy(v)
			end
		end
	end

	printstyled(
		"\e[2K\e[1GDone! Data saved at " * FilePathOut * "\n", color=:green
	)

end

@doc raw"""
...
"""
function RunHFRecord(
	Phase::String,						# Mean field phase
	t::Float64,							# Hopping amplitude
	U::Float64,							# Local repulsion
	V::Float64,							# Non-local attraction
	Lx::Int64,							# Lattice size
	δ::Float64,							# Doping
	β::Float64,							# Inverse temperature
	p::Int64,							# Number of iterations
	Δv::Dict{String,Float64},			# Tolerance on magnetization
	Δn::Float64,						# Tolerance on density
	gg::Vector{Float64};				# Mixing parameter
	Syms::Vector{String}=["d"],			# Gap function symmetries
	DirPathOut::String="",				# Output file
	RenormalizeBands::Bool=true		# Conditional renormalization of t
)::Dict{Float64,Dict{String,Vector{Float64}}}

	L = [Lx,Lx]

	# Get Hartree Fock Parameters labels
	HFPs = GetHFPs(Phase)

	# Initialize model parameters
	Parameters::Dict{String,Float64} = Dict([
		"t" => t,
		"U" => U,
		"V" => V
	])

	printstyled(
		"Recording HF at t=$t\t U=$U\t V=$V\t L=$L\t β=$β\t δ=$δ\n", 
		color=:yellow
	)

	Record::Dict{Float64,Dict{String,Vector{Float64}}} = Dict([])
	for g in gg

		# Record routine
		HFResults = RunHFAlgorithm(
			Phase,Parameters,L,0.5+δ,β,
			p,Δv,Δn,g;
			verbose=true,
			record=true,
			RenormalizeBands
		)
		ΔT::Float64 = HFResults[3]
		gRecord::Dict{String,Vector{Float64}} = Dict([
			key => HFResults[5][key] for key in HFPs
		])

		# Write record on matrix
		RecordMatrix::Matrix{Float64} = zeros(
			length(values(
				gRecord[ HFPs[1] ]
			)),
			length( keys(HFPs) )
		)
		for (k,key) in enumerate(HFPs)
			RecordMatrix[:,k] = gRecord[key]
		end

		# Write on file
		if DirPathOut != ""
			FilePathOut = DirPathOut * "g=$(g).txt"

			# File initialization
			Header = "# $(HFPs) [calculated @ $(now())]\n"
			write(FilePathOut, Header)

			# Append recorded matrix
			open(FilePathOut, "a") do io
				writedlm(io, RecordMatrix, ';')
			end
			printstyled(
				"\e[2K\e[1GDone! Data saved at " * FilePathOut *
				"\n", color=:green
			)
		end
		Record[g] = gRecord
	end
	return Record
end

# Main run
function main()
	DirPathOut = PROJECT_ROOT * "/simulations/" *
		Mode * "/Setup=$(Setup)/Phase="
	if in(Mode, ["scan", "heatmap"])
		if in(Phase, ["SU-Singlet", "FakeSU-Singlet", "SU-Triplet", "FakeSU-Triplet"])
			FilePathOut = DirPathOut * Phase * "/Syms=$(Syms...).txt"
		else
			FilePathOut = DirPathOut * Phase * ".txt"
		end
		mkpath(dirname(FilePathOut))
		RunHFScan(
			Phase,
			tt,UU,VV,
			LL,δδ,ββ,
			p,Δv,Δn,g;
			Syms,
			FilePathOut,
			RenormalizeBands
		)
	elseif Mode=="record-g"
		RunHFRecord(
			Phase,
			t,U,V,
			L,δ,β,
			p,Δv,Δn,gg;
			Syms,
			DirPathOut,
			RenormalizeBands
		)
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
