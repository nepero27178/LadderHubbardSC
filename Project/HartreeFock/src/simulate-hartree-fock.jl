#!/usr/bin/julia
PROJECT_ROOT = @__DIR__
include(PROJECT_ROOT * "/setup/simulations-setup.jl")
include(PROJECT_ROOT * "/modules/reciprocal-space-algorithm.jl")

using DelimitedFiles
using Dates

function RunHFRoutine(
    t::Float64,
    LL::Vector{Int64},
    δδ::Vector{Float64},
    ββ::Vector{Float64},
    p::Int64,
    Δm::Float64,
    Δn::Float64,
    FilePathOut::String
)

    Header = "# Lx, β, δ, m, Q, ΔT [calculated @ $(now())]\n"
    write(FilePathOut, Header)

    ResultsTable = zeros( length(LL) * length(ββ) * length(δδ), 6 )
    i = 1    
    for Lx in LL, β in ββ, δ in δδ
        L = [Lx, Lx]
        printstyled("\e[2K\e[1GRunning HF for Lx=$Lx, β=$β, δ=$δ", color=:yellow)
        ResultsTable[i,1:3] .= [Lx,β,δ]
        ResultsTable[i,4:6] .= RunHFAlgorithm(t,L,0.5+δ,β,p,Δm,Δn)
        i += 1
    end
    
    open(FilePathOut, "a") do io
        writedlm(io, ResultsTable, ',')
    end
    printstyled("\e[2K\e[1GDone! Data saved at $(FilePathOut)\n", color=:green)    

end

function main()
    DirPathOut = PROJECT_ROOT * "/../simulations/hubbard/p=$(p)_Δm=$(Δm)_Δn=$(Δn)/"
    mkpath(DirPathOut)
    FilePathOut = DirPathOut * "t=$(t)_ββ=$(ββ).txt"
    RunHFRoutine(t,LL,δδ,ββ,p,Δm,Δn,FilePathOut)
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
