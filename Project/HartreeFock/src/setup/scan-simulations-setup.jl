#!/usr/bin/julia
SetupFilePath = @__FILE__

# All possible simulations
AllPhases = ["AF", "SU/Singlet", "SU/Triplet"]
AllSingletSyms = ["s", "s*", "d"]
AllTripletSyms = ["px", "py", "p+", "p-"]

Phase = "AF"    # Choose your phase
Model = "Renormalized-AF"
Setup = "B"  # Choose your setup

if !in(Setup, ["Test", "A", "B"])
    @error "Invalid setup, please modify at: " * SetupFilePath
    exit()
elseif Setup=="Test"
    # TEST-SETUP
    tt = [1.0]
    UU = [U for U in 0.1:0.1:10.0]
    VV = [0.0]
    LL = [16]
    δδ = [δ for δ in 0.0:0.1:0.45]
    ββ = [100.0, Inf]
    p = 100
	Δv::Dict{String,Float64} = Dict([
	    "m" => 1e-3,
	    "w0" => 1e-3,
	    "wp" => 1e-3,
	])
    Δn = 1e-2
    g = 0.5
elseif Setup=="A"
    tt = [1.0]
    UU = [4.0,10.0,20.0]
    VV = [V for V in 0.0:0.1:2.0]
    LL = [32]
    δδ = [δ for δ in 0.0:0.05:0.45]
    ββ = [100.0, Inf]
    p = 100
	Δv::Dict{String,Float64} = Dict([
	    "m" => 1e-3,
	    "w0" => 1e-3,
	    "wp" => 1e-3
	])
    Δn = 1e-2
    g = 0.5
elseif Setup=="B"
    tt = [1.0]
    UU = [U for U in 4.0:4.0:20.0]
    VV = [V for V in 0.0:0.2:4.0]
    LL = [32,48,64]
    δδ = [δ for δ in 0.0:0.05:0.45]
    ββ = [10.0, 50.0, 100.0, Inf]
    p = 1000
	Δv::Dict{String,Float64} = Dict([
	    "m" => 1e-4,
	    "w0" => 1e-4,
	    "wp" => 1e-4
	])
    Δn = 1e-2
    g = 0.49
end
