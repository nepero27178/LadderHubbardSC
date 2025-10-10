#!/usr/bin/julia
SetupFilePath = @__FILE__

# All possible simulations
AllPhases = ["AF", "SU/Singlet", "SU/Triplet"]
AllSingletSyms = ["s", "s*", "d"]
AllTripletSyms = ["px", "py", "p+", "p-"]

Phase = "AF"    # Choose your phase
Setup = "Test"  # Choose your setup

if !in(Setup, ["Test"])
    @error "Invalid setup, please modify at: " * SetupFilePath
    exit()

elseif Setup=="Test"
    # TEST-SETUP
    tt = [1.0]
    UU = [10.0]
    VV = [V for V in 0.8:0.2:1.0]
    LL = [16]
    δδ = [δ for δ in 0.0:0.1:0.2]
    ββ = [100.0, Inf]
    p = 100
	Δv::Dict{String,Float64} = Dict([
	    "m" => 1e-1,
	    "w^0" => 1e-1,
	    "w^pi" => 1e-1
	])
    Δn = 1e-2
    g = 0.5
end
