#!/usr/bin/julia
PROJECT_ROOT = @__DIR__
include(PROJECT_ROOT * "/setup/simulations-setup.jl")
include(PROJECT_ROOT * "/modules/reciprocal-space-algorithm.jl")

using DelimitedFiles
using Dates

function RunHFRoutine(
    UU::Vector{Float64},
    LL::Vector{Int64},
    δδ::Vector{Float64},
    ββ::Vector{Float64},
    p::Int64,
    Δm::Float64,
    Δn::Float64,
    FilePathOut::String
)
    
    tt = 1 ./ UU   # Convert to hopping variable

    Header = "# U, t, Lx, β, δ, m, Q, ΔT [calculated @ $(now())]\n"
    write(FilePathOut, Header)

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
        ResultsTable[i,6:8] .= RunHFAlgorithm(t,L,0.5+δ,β,p,Δm,Δn)
        i += 1
    end
    
    open(FilePathOut, "a") do io
        writedlm(io, ResultsTable, ',')
    end
    printstyled("\e[2K\e[1GDone! Data saved at $(FilePathOut)\n", color=:green)    

end

function main()
    DirPathOut = PROJECT_ROOT *
        "/../simulations/hubbard/p=$(p)_Δm=$(Δm)_Δn=$(Δn)/"
    mkpath(DirPathOut)
    FilePathOut = DirPathOut * "Setup=$(Setup).txt"
    RunHFRoutine(UU,LL,δδ,ββ,p,Δm,Δn,FilePathOut)
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
