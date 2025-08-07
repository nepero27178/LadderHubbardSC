#!/usr/bin/julia
SetupFilePath = @__FILE__

Setup="Test" # Choose your setup

if !in(Setup, ["Test", "A", "B"])
    @error "Invalid setup, please modify at: " * SetupFilePath
    exit()

elseif Setup=="Test"

    # TEST-SETUP
    UU = [U for U in 0.5:0.5:10.0]
    LL = [32] # [2^x for x in 4:7]
    δδ = [δ for δ in 0.0:0.1:0.3]
    ββ = [Inf] # [β for β in 10.0:20.0:50.0]
    p = 100
    Δm = 1e-4
    Δn = 1e-2

elseif Setup=="A"
    # SETUP A ...
    @error "Empty setup selected, please modify at: " * SetupFilePath

elseif Setup=="B"
    # SETUB B ...
    @error "Empty setup selected, please modify at: " * SetupFilePath

end
