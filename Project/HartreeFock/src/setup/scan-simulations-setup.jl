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
Setup = "B[256]"  # Choose your setup #TODO Use readline()
AvailableSetups = ["Test256", "A[128,256]", "A[128,256]-PureHopping", "B[256]", "B[256]-PureHopping"]

if !in(Setup, AvailableSetups)
    @error "Invalid setup, please modify at: " * SetupFilePath
    exit()
elseif Setup=="Test256"
    # TEST-SETUP
    tt = [1.0]
    UU = [U for U in 0.0:1.0:20.0]
    VV = [0.0]
    LL = [256]
    δδ = [δ for δ in 0.0:0.05:0.45]
    ββ = [100.0, Inf]
    p = 100
	Δv::Dict{String,Float64} = Dict([
	    "m" => 1e-4,
	    "w0" => 1e-4,
	    "wp" => 1e-4,
	])
    Δn = 1e-2
    g = 0.5
elseif in(Setup, ["A[128,256]", "A[128,256]-PureHopping"])
    tt = [1.0]
    UU = [4.0,10.0,20.0]
    VV = [V for V in 0.0:0.1:3.0]
    LL = [128,256]
    δδ = [δ for δ in 0.0:0.05:0.45]
    ββ = [50.0, 100.0]
    p = 100
	Δv::Dict{String,Float64} = Dict([
	    "m" => 1e-4,
	    "w0" => 1e-4,
	    "wp" => 1e-4
	])
    Δn = 1e-2
    g = 0.5
    RenormalizeHopping::Bool = true
    if Setup=="A[128,256]-PureHopping"
        RenormalizeHopping::Bool = false
    end
elseif in(Setup, ["B[256]", "B[256]-PureHopping"])
    tt = [1.0]
    UU = [4.0,12.0,20.0]
    VV = [V for V in 0.0:0.1:3.0]
    LL = [256]
    δδ = [δ for δ in 0.0:0.05:0.45]
    ββ = [100.0]
    p = 100
	Δv::Dict{String,Float64} = Dict([
	    "m" => 1e-4,
	    "w0" => 1e-4,
	    "wp" => 1e-4
	])
    Δn = 1e-2
    g = 0.5
    RenormalizeHopping::Bool = true
    if Setup=="B[256]-PureHopping"
        RenormalizeHopping::Bool = false
    end
end
