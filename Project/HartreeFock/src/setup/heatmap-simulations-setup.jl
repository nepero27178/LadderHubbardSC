#!/usr/bin/julia
SetupFilePath = @__FILE__

# All possible simulations
AllPhases = [
    "AF",               # Renormalized AntiFerromagnet
    "AF*",              # As for AF, but with pure Hopping
    "SU/Singlet",       # Singlet superconductor
    "SU/Triplet"        # Triplet superconductor
]
AllSingletSyms = ["s", "s*", "d"]
AllTripletSyms = ["px", "py", "p+", "p-"]

Phase = "AF*"    # Choose your phase
if !in(Phase, AllPhases)
    @error "Invalid phase, please modify at: " * SetupFilePath
    exit()
end
# Model = "Renormalized-AF"
Setup = "Test[50]"  # Choose your setup #TODO Use readline()
AvailableSetups = [
    "Test[50]",
    "A[128]",
    "B[128]",
    "C[128]"
]

RenormalizeHopping::Bool = true
if Phase=="AF*"
    RenormalizeHopping::Bool = false
end

if !in(Setup, AvailableSetups)
    @error "Invalid setup, please modify at: " * SetupFilePath
    exit()
elseif Setup=="Test[50]"
    # TEST-SETUP
    tt = [1.0]
    UU = [0.0, 10.0, 20.0]
    VV = [1.0, 2.0, 3.0]
    LL = [80]
    δδ = [0.1, 0.2, 0.3]
    ββ = [100.0]
    p = 100
	Δv::Dict{String,Float64} = Dict([
	    "m" => 1e-2,
	    "w0" => 1e-2,
	    "wp" => 1e-2
	])
    Δn = 1e-2
    g = 0.5
elseif Setup=="A[128]"
    tt = [1.0]
    UU = [U for U in 0.0:1.0:20.0]
    VV = [V for V in 0.0:0.1:4.0]
    LL = [128]
    δδ = [0.2]
    ββ = [100.0]
    p = 100
	Δv::Dict{String,Float64} = Dict([
	    "m" => 1e-4,
	    "w0" => 1e-4,
	    "wp" => 1e-4
	])
    Δn = 1e-2
    g = 0.5
elseif Setup=="B[128]"
    tt = [1.0]
    UU = [4.0]
    VV = [V for V in 0.0:0.1:4.0]
    LL = [128]
    δδ = [δ for δ in 0.0:0.01:0.49]
    ββ = [100.0]
    p = 100
	Δv::Dict{String,Float64} = Dict([
	    "m" => 1e-4,
	    "w0" => 1e-4,
	    "wp" => 1e-4
	])
    Δn = 1e-2
    g = 0.5
elseif Setup=="C[128]"
    tt = [1.0]
    UU = [12.0]
    VV = [V for V in 0.0:0.1:4.0]
    LL = [128]
    δδ = [δ for δ in 0.0:0.01:0.49]
    ββ = [100.0]
    p = 100
	Δv::Dict{String,Float64} = Dict([
	    "m" => 1e-4,
	    "w0" => 1e-4,
	    "wp" => 1e-4
	])
    Δn = 1e-2
    g = 0.5
end
