#!/usr/bin/julia
PROJECT_ROOT = @__DIR__
PROJECT_ROOT *= "/.."   # Up to the effective root
include(PROJECT_ROOT * "/src/setup/simulations-setup.jl")
include(PROJECT_ROOT * "/src/modules/reciprocal-space-algorithm.jl")

using DelimitedFiles
using Dates

@doc raw"""
function RunHFRoutine(
    UU::Vector{Float64},
    LL::Vector{Int64},
    δδ::Vector{Float64},
    ββ::Vector{Float64},
    p::Int64,
    Δm::Float64,
    Δn::Float64,
    g::Float64;
    FilePathOut::String=""
)

Returns: none if `FilePathOut` is specified, `ResultsTable::Matrix{Float64}` if
`FilePathOut` is unspecified.

`RunHFRoutine` takes as input `UU` (vector of local repulsions), `LL` (vector of
the square lattice dimensions), `δδ` (vector of dopings with respect to the half
filled lattice), `ββ` (vector of inverse temperatures), `p` (maximum number of
HF iterations), `Δm` (tolerance on magnetization) and `Δn` (tolerance on density
in chemical potential estimation), `g` (mixing parameter). It performs an 
iterative HF analysis over a  sequence of 2D square lattices for all the 
possible combinations of the specified parameters. Check the source files in the
`/modules` folder for more informations on the algorithm.
"""
function RunHFRoutine(
    UU::Vector{Float64},        # Local repulsion
    LL::Vector{Int64},          # Lattice size
    δδ::Vector{Float64},        # Doping
    ββ::Vector{Float64},        # Inverse temperature
    p::Int64,                   # Number of iterations
    Δm::Float64,                # Tolerance on magnetization
    Δn::Float64,                # Tolerance on density
    g::Float64;                 # Mixing parameter
    FilePathOut::String=""      # Output file
)
    
    tt = 1 ./ UU   # Convert to hopping variable

    if FilePathOut != ""
       Header = "# U, t, Lx, β, δ, m, Q, ΔT [calculated @ $(now())]\n"
        write(FilePathOut, Header)
    end

    ResultsTable = zeros( length(UU) * length(LL) * length(ββ) * length(δδ), 8 )
    i = 1    
    for (u,U) in enumerate(UU), Lx in LL, β in ββ, δ in δδ
        L = [Lx, Lx]
        t = tt[u]
        printstyled(
            "\e[2K\e[1GRunning HF for U=$U, Lx=$Lx, β=$β, δ=$δ", 
            color=:yellow
        )
        ResultsTable[i,1:5] .= [U,t,Lx,β,δ]
        ResultsTable[i,6:8] .= RunHFAlgorithm(t,L,0.5+δ,β,p,Δm,Δn,g)
        i += 1
    end
    
    if FilePathOut != ""
        open(FilePathOut, "a") do io
            writedlm(io, ResultsTable, ',')
        end
        printstyled(
            "\e[2K\e[1GDone! Data saved at " * FilePathOut * "\n",
            color=:green
        )    
    elseif FilePathOut==""
        printstyled(
            "\e[2K\e[1GDone!\n",
            color=:green
        )
        return ResultsTable
    end

end

# Main run
function main()
    DirPathOut = PROJECT_ROOT * "/simulations/hubbard/p=$(p)_Δm=$(Δm)_Δn=$(Δn)/"
    mkpath(DirPathOut)
    FilePathOut = DirPathOut * "Setup=$(Setup).txt"
    RunHFRoutine(UU,LL,δδ,ββ,p,Δm,Δn,g;FilePathOut)
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
