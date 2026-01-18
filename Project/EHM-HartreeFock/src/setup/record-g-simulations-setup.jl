#!/usr/bin/julia
SetupFilePath = @__FILE__

# All possible simulations
AllPhases = [
	"AF",               # Renormalized AntiFerromagnet
	"FakeAF",              # As for AF, but with pure Hopping
	"SU-Singlet",       # Singlet superconductor
	"SU-Triplet"        # Triplet superconductor
]
AllSingletSyms = ["s", "S", "d"]
AllTripletSyms = ["px", "py", "p+", "p-"]

Phase = "AF"    # Choose your phase
if !in(Phase, AllPhases)
	@error "Invalid phase, please modify at: " * SetupFilePath
	exit()
end
# Model = "Renormalized-AF"
Setup = "Test16"  # Choose your setup #TODO Use readline()
AvailableSetups = [
	"Test16",
	"A32",
	"A256",
	"B128",
	"B256"
]

RenormalizeBands::Bool = true
if Phase=="FakeAF"
	RenormalizeBands::Bool = false
end

if !in(Setup, AvailableSetups)
	@error "Invalid setup, please modify at: " * SetupFilePath
elseif Setup=="Test16"
	# A point, L=32
	t = 1.0
	U = 10.0
	V = 2.0
	L = 16
	δ = 0.4
	β = 100.0
	p = 20
	Δv::Dict{String,Float64} = Dict([
		"m" => 1e-2,
		"w0" => 1e-2,
		"wp" => 1e-2,
	])
	Δn = 1e-2
	gg = [0.05, 0.5, 0.95]
elseif Setup=="A32"
	# A point, L=32
	t = 1.0
	U = 10.0
	V = 5.3
	L = 32
	δ = 0.4
	β = 100.0
	p = 1000
	Δv::Dict{String,Float64} = Dict([
		"m" => 1e-4,
		"w0" => 1e-3,
		"wp" => 1e-4,
	])
	Δn = 1e-2
	gg = [g for g in 0.04:0.04:0.95]
elseif Setup=="A128"
	# A point, L=128
	t = 1.0
	U = 10.0
	V = 5.3
	L = 128
	δ = 0.4
	β = 100.0
	p = 1000
	Δv::Dict{String,Float64} = Dict([
		"m" => 1e-4,
		"w0" => 1e-3,
		"wp" => 1e-4,
	])
	Δn = 1e-2
	gg = [g for g in 0.04:0.04:0.95]
elseif Setup=="B32"
	# B point, L=32
	t = 1.0
	U = 4.0
	V = 0.8
	L = 32
	δ = 0.4
	β = Inf
	p = 1000
	Δv::Dict{String,Float64} = Dict([
		"m" => 1e-4,
		"w0" => 1e-3,
		"wp" => 1e-4,
	])
	Δn = 1e-2
	gg = [g for g in 0.04:0.04:0.95]
elseif Setup=="B128"
	# B point, L=128
	t = 1.0
	U = 4.0
	V = 0.8
	L = 128
	δ = 0.4
	β = Inf
	p = 1000
	Δv::Dict{String,Float64} = Dict([
		"m" => 1e-4,
		"w0" => 1e-3,
		"wp" => 1e-4,
	])
	Δn = 1e-2
	gg = [g for g in 0.04:0.04:0.95]
end
