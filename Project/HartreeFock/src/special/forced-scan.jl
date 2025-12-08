#!/usr/bin/julia
using DelimitedFiles
using DataFrames

# Includer
PROJECT_ROOT = @__DIR__
PROJECT_ROOT *= "/../.." # Up to the effective root
include(PROJECT_ROOT * "/src/modules/methods-simulating.jl")

function RunHFScan_AF_ForceHoppingShift(
	DF::DataFrame,						# Input DataFrame to be reproduced
	FilePathOut::String="",				# Output file
)

	Phase = "AF"
	p = 100
	Δv = Dict([
		"m" => 1e-4,
		"w0" => 1e-4,
		"wp" => 1e-4
	])
	Δn = 1e-2
	g = 0.5

	# Get Hartree Fock Parameters labels
    HFPs = GetHFPs(Phase)

	# File coditional initialization (otherwise, just append)
	if FilePathOut != ""
		Header = "t;U;V;Lx;β;δ;v;Q;ΔT;μ\n"
		write(FilePathOut, Header)
	end

	# HF iterations
	i = 1
	Iterations = size(DF,1)
	for Point in eachrow(DF)

		# Initializer
		v0::Dict{String,Float64} = eval(Meta.parse(Point.v))

		# Model parameters
		w0 = eval(Meta.parse(Point.v))["w0"]
		t = Point.t - w0 * Point.V
		U = Point.U
		V = Point.V

		Parameters::Dict{String,Float64} = Dict([
			"t" => t,
			"U" => U,
			"V" => V
		])

		Lx = Int64(Point.Lx)
		β = Point.β
		δ = Point.δ
		L = [Lx, Lx]
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
			v0i=v0,
			RenormalizeHopping=false
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
function ForcedScan()
	Setup = "B[256]"
	FilePathIn = PROJECT_ROOT * "/simulations/scan/Setup=$(Setup)/AF.txt"
	DataIn = open(FilePathIn) do io
		readdlm(FilePathIn, ';', comments=false, '\n')
	end
	DF = DataFrame(DataIn[2:end,:], DataIn[1,:])

	FilePathOut = replace(FilePathIn, "scan" => "special/forced-scan")
	mkpath(dirname(FilePathOut))

	RunHFScan_AF_ForceHoppingShift(DF,FilePathOut)
end

if abspath(PROGRAM_FILE) == @__FILE__
	ForcedScan()
end