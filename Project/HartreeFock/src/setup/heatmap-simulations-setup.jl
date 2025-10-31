#!/usr/bin/julia
SetupFilePath = @__FILE__

# All possible simulations
AllPhases = ["AF", "SU/Singlet", "SU/Triplet"]
AllSingletSyms = ["s", "s*", "d"]
AllTripletSyms = ["px", "py", "p+", "p-"]

Phase = "AF"    # Choose your phase
if !in(Phase, AllPhases)
    @error "Invalid phase, please modify at: " * SetupFilePath
    exit()
end
Model = "Renormalized-AF"
Setup = "B128"  # Choose your setup
RenormalizeHopping::Bool = false

if !in(Setup, ["Test", "A128", "B128"])
    @error "Invalid setup, please modify at: " * SetupFilePath
    exit()
elseif Setup=="Test"
    # TEST-SETUP
    tt = [1.0]
    UU = [U for U in 0.0:1.0:3.0]
    VV = [V for V in 0.0:10.0:20.0]
    LL = [64]
    δδ = [0.2]
    ββ = [100.0]
    p = 100
	Δv::Dict{String,Float64} = Dict([
	    "m" => 1e-2,
	    "w0" => 1e-2,
	    "wp" => 1e-2
	])
    Δn = 1e-2
    g = 0.5
elseif Setup=="A128"
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
elseif Setup=="B128"
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
end
