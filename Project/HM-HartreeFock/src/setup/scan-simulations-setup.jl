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
Setup = "A[128]"  # Choose your setup #TODO Use readline()
AvailableSetups = [
	"Test[80]",
	"A[128]"
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
elseif Setup=="A[128]"
    #...
end
