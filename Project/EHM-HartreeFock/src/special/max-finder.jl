#!/usr/bin/julia
using DelimitedFiles
using Dates

# Includer
PROJECT_ROOT = @__DIR__
PROJECT_ROOT *= "/../.." # Up to the effective root
include(PROJECT_ROOT * "/src/modules/methods-simulating.jl")

function RunHFScan_AF_MaxFinder(
	t::Float64,							# Hopping amplitude
	U::Float64,							# Local repulsion
	Lx::Int64,							# Lattice size
	β::Float64,							# Inverse temperature
	VV::Vector{Float64},					# Non-local attraction
	uu::Vector{Float64},					# Lower boundary
	ll::Vector{Float64},					# Upper boudary
	Steps::Int64,						# Number of steps
	p::Int64,							# Number of iterations
	Δv::Dict{String,Float64},			# Tolerance on magnetization
	Δn::Float64,							# Tolerance on density
	g::Float64;							# Mixing parameter
	Phase::String="AF",					# Model phase
	FilePathOut::String="",				# Output file
	RenormalizeHopping::Bool=true		# Conditional renormalization of t
)
	@warn "This function is going to be deprecated."

	# Warn user of memory-heavy simulations detection
	Iterations = length(VV) * Steps
	if FilePathOut=="" && Iterations > 200
		@warn "No output file specified and more than 200 simulations " *
			"request detected. Simulations results are going to be stored " *
			"in your memory. Consider specifying a `FilePathOut` and " *
			"unloading your memory."
	end

	# Get Hartree Fock Parameters labels
	HFPs = GetHFPs("AF")

	# File coditional initialization (otherwise, just append)
	if FilePathOut != ""
			Header = "t;U;V;Lx;β;δ;v;Q;ΔT\n"
		write(FilePathOut, Header)
	end

	# HF iterations
	i = 1
	for (j,V) in enumerate(VV),

		δδ = range(ll[j], uu[j], Steps)
		Parameters::Dict{String,Float64} = Dict([
			"t" => t,
			"U" => U,
			"V" => V
		])

		L = [Lx, Lx]

		for δ in δδ
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
				RenormalizeHopping
			)

			v::Dict{String,Float64} = Dict([
				key => HFResults[1][key] for key in HFPs
			])
			Qs::Dict{String,Float64} = Dict([
				key => HFResults[2][key] for key in HFPs
			])
			ResultsVector = hcat(ResultsVector[:,3:end], [v Qs HFResults[3]])

			i += 1

			# Append to initialized or existing file
			if FilePathOut != ""
				open(FilePathOut, "a") do io
					writedlm(io, ResultsVector, ';')
				end
			end
		end
	end

	printstyled(
		"\e[2K\e[1GDone! Data saved at " * FilePathOut * "\n", color=:green
	)
end

# Main run
function MaxFinder()
	DirPathOut = PROJECT_ROOT * "/simulations/special/max-finder/"
	FilePathOut = DirPathOut * "$(today())-MaxFinder-AF.txt"
	mkpath(dirname(DirPathOut))

	# Model setup
	t::Float64 = 1.0
	U::Float64 = 4.0
	β::Float64 = 100.0
	Lx::Int64 = 256

	# Algorithm setup
	p::Int64 = 100
	Δv::Dict{String,Float64} = Dict([
		"m" => 1e-4,
		"w0" => 1e-4,
		"wp" => 1e-4
	])
	Δn::Float64 = 1e-2
	g::Float64 = 0.5

	# Special selection
	VV::Vector{Float64} = [V for V in 0.0:0.04:4.0]

	# Rough lab fit results
	a::Float64 = -0.05
	b::Float64 = 0.25
	n::Float64 = 0.40

	# Upper and lower boundaries
	uu::Vector{Float64} = a .+ b .* VV.^ n .+ 0.03
	ll::Vector{Float64} = a .+ b .* VV.^ n .- 0.03
	Steps::Int64 = 30 # Arbitrary

	# Run
	RunHFScan_AF_MaxFinder(
		t,U,Lx,β,
		VV,uu,ll,Steps,
		p,Δv,Δn,g;
		FilePathOut
	)
end

if abspath(PROGRAM_FILE) == @__FILE__
	SpecialScan()
end
