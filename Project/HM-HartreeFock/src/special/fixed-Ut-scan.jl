#!/usr/bin/julia
using DelimitedFiles
using DataFrames

# Includer
PROJECT_ROOT = @__DIR__
PROJECT_ROOT *= "/../.." # Up to the effective root
include(PROJECT_ROOT * "/src/setup/graphic-setup.jl")
include(PROJECT_ROOT * "/src/modules/methods-simulating.jl")
include(PROJECT_ROOT * "/src/modules/methods-plotting.jl")

function RunHFScan_AF_FixedUt(
	Phase::String,						# Mean field phase
	xx::Vector{Float64},				# Ratios
	U::Float64,							# Fized local repulsion
	LL::Vector{Int64},					# Lattice size
	δδ::Vector{Float64},				# Doping
	ββ::Vector{Float64},				# Inverse temperature
	p::Int64,							# Number of iterations
	Δv::Dict{String,Float64},			# Tolerance on magnetization
	Δn::Float64,						# Tolerance on density
	g::Float64;							# Mixing parameter
	FilePathOut::String="",				# Output file
	InitializeFile::Bool=true,			# Initialize file at FilePathOut
)

	# Warn user of memory-heavy simulations detection
	Iterations = length(xx) * length(LL) * length(ββ) * length(δδ)
	if FilePathOut == "" && Iterations > 200
		@warn "No output file specified and more than 200 simulations " *
			"request detected. Simulations results are going to be stored " *
			"in your memory. Consider specifying a `FilePathOut` and " *
			"unloading your memory."
	end

	# Get Hartree Fock Parameters labels
	HFPs = GetHFPs(Phase)

	# File coditional initialization (otherwise, just append)
	if FilePathOut != "" && InitializeFile
		Header = "x;t;U;Lx;β;δ;v;Q;ΔT;μ\n"
		write(FilePathOut, Header)
	end

	# Initializers
	v0::Dict{String,Float64} = Dict([
		key => 0.5 for key in HFPs
	])

	# HF iterations
	i = 1
	for x in xx,
		Lx in LL,
		β in ββ,
		δ in δδ

		t = U/x
		Parameters::Dict{String,Float64} = Dict([
			"t" => t,
			"U" => U
		])

		L = [Lx, Lx]
		printstyled(
			"\e[2K\e[1GRun ($i/$Iterations): " *
			"$Phase HF at x=$x, t=$t, U=$U, L=$Lx, β=$β, δ=$δ", 
			color=:yellow
		)
		ResultsVector::Matrix{Any} = [0 0]  # Dummy placeholder
		ResultsVector = hcat(ResultsVector, [x t U Lx β δ])

		# Run routine, all positional arguments here must be false
		HFResults = RunHFAlgorithm(
			Phase,Parameters,L,0.5+δ,β,
			p,Δv,Δn,g;
			#v0i=v0
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

# Main run
function main()

	xx = [x for x in 0.5:0.5:20.0]
	U = 10.0
	LL = [128]
	δδ = [δ for δ in 0.0:0.05:0.45]
	ββ = [100.0]
	p = 100
	Δv::Dict{String,Float64} = Dict([
		"m" => 1e-4,
	])
	Δn = 1e-2
	g = 0.5

	FilePathOut = PROJECT_ROOT * "/simulations/special/fixed-Ut-scan/vary-U/AF.txt"
	mkpath(dirname(FilePathOut))
	RunHFScan_AF_FixedUt(
		"AF",
		xx,U,
		LL,δδ,ββ,
		p,Δv,Δn,g;
		FilePathOut
	)
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
