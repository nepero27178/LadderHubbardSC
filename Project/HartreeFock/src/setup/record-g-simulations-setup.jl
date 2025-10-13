#!/usr/bin/julia
SetupFilePath = @__FILE__

# All possible simulations
AllPhases = ["AF", "SU/Singlet", "SU/Triplet"]
AllSingletSyms = ["s", "s*", "d"]
AllTripletSyms = ["px", "py", "p+", "p-"]

Phase = "AF"    # Choose your phase
Model = "Renormalized-AF"
Setup = "A"  # Choose your setup

if !in(Setup, ["A"])
    @error "Invalid setup, please modify at: " * SetupFilePath
    exit()
elseif Setup=="A"
    # TEST-SETUP
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
end
