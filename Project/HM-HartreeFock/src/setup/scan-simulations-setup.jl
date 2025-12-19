#!/usr/bin/julia
SetupFilePath = @__FILE__

# All possible simulations
AllPhases = [
	"AF",               # Renormalized AntiFerromagnet
]

Phase = "AF"    # Choose your phase
if !in(Phase, AllPhases)
	@error "Invalid phase, please modify at: " * SetupFilePath
	exit()
end
# Model = "Renormalized-AF"
Setup = "A[128]_t=1.0"  # Choose your setup #TODO Use readline()
AvailableSetups = [
	"Test[80]",
	"A[128]_t=1.0",
	"A[128]_t=0.5",
	"A[128]_t=2.0"
]

if !in(Setup, AvailableSetups)
	@error "Invalid setup, please modify at: " * SetupFilePath
elseif Setup=="Test[80]"
	# TEST-SETUP
	tt = [1.0]
	UU = [10.0]
	LL = [80]
	δδ = [0.1, 0.2, 0.3]
	ββ = [100.0]
	p = 20
	Δv::Dict{String,Float64} = Dict([
		"m" => 1e-3,
	])
	Δn = 1e-2
	g = 0.5
elseif Setup=="A[128]_t=1.0"
	tt = [1.0]
	UU = [U for U in 1.0:1.0:20.0]
	LL = [128]
	δδ = [δ for δ in 0.0:0.05:0.45]
	ββ = [100.0]
	p = 100
	Δv::Dict{String,Float64} = Dict([
		"m" => 1e-4,
	])
	Δn = 1e-2
	g = 0.5
elseif Setup=="A[128]_t=0.5"
	tt = [0.5]
	UU = [U for U in 0.5:0.5:10.0]
	LL = [128]
	δδ = [δ for δ in 0.0:0.05:0.45]
	ββ = [100.0]
	p = 100
	Δv::Dict{String,Float64} = Dict([
		"m" => 1e-4,
	])
	Δn = 1e-2
	g = 0.5
elseif Setup=="A[128]_t=2.0"
	tt = [2.0]
	UU = [U for U in 2.0:2.0:40.0]
	LL = [128]
	δδ = [δ for δ in 0.0:0.05:0.45]
	ββ = [100.0]
	p = 100
	Δv::Dict{String,Float64} = Dict([
		"m" => 1e-4,
	])
	Δn = 1e-2
	g = 0.5
end
