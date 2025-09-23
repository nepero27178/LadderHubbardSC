#!/usr/bin/julia
SetupFilePath = @__FILE__

AllSyms = ["s", "s*", "px", "py", "d"]  # Gap function possible symmetries

HMSyms = ["d"]
HMSymsStr = ""
for Sym in HMSyms
    global HMSymsStr *= "$(Sym)-"
end
# Pop last "-" character
global HMSymsStr = HMSymsStr[1:end-1]

Setup="C" # Choose your setup

if !in(Setup, ["Test", "A", "B", "C"])
    @error "Invalid setup, please modify at: " * SetupFilePath
    exit()

elseif Setup=="Test"

    # TEST-SETUP
    UU = [10.0]
    VV = [V for V in 0.8:0.2:1.0]
    LL = [16]
    δδ = [δ for δ in 0.0:0.1:0.2]
    ββ = [100.0, Inf]
    p = 100
    Δm = Dict([
        Sym => 1e-1 for Sym in AllSyms
    ])
    Δn = 1e-2
    g = 0.5

elseif Setup=="A"
    # SETUP A ...
    UU = [4.0, 10.0]                    # Local repulsions
    VV = [V for V in 0.1:0.1:1.0]		# Non-local attractions
    LL = [2^5]                          # Lattice sizes
    δδ = [δ for δ in 0.0:0.05:0.45]     # Dopings
    ββ = [10.0, Inf]				    # Inverse temperatures
    p = 100                             # Max number of iterations
    Δm = Dict([                         # Tolerance on each symmetry
        Sym => 1e-3 for Sym in AllSyms
    ])
    Δn = 1e-2                           # Tolerance on density
    g = 0.5                             # Mixing parameter

elseif Setup=="B"
    # SETUB B ...
    # SETUP A ...
    UU = [10.0]		                    # Local repulsions
    VV = [V for V in 0.0:0.25:2.5]		# Non-local attractions
    LL = [2^5]                          # Lattice sizes
    δδ = [δ for δ in 0.0:0.05:0.45]     # Dopings
    ββ = [Inf]						    # Inverse temperatures
    p = 100                             # Max number of iterations
    Δm = Dict([                         # Tolerance on each symmetry
        Sym => 1e-3 for Sym in AllSyms
    ])
    Δn = 1e-2                           # Tolerance on density
    g = 0.5                             # Mixing parameter

elseif Setup=="C"
    UU = [10.0]		                    # Local repulsions
    VV = [V for V in 1.0:0.1:2.0]		# Non-local attractions
    LL = [2^5]                          # Lattice sizes
    δδ = [δ for δ in 0.3:0.01:0.45]     # Dopings
    ββ = [Inf]						    # Inverse temperatures
    p = 100                             # Max number of iterations
    Δm = Dict([                         # Tolerance on each symmetry
        Sym => 1e-3 for Sym in AllSyms
    ])
    Δn = 1e-2                           # Tolerance on density
    g = 0.5                             # Mixing parameter
end
