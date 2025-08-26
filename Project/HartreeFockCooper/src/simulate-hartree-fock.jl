#!/usr/bin/julia
PROJECT_ROOT = @__DIR__
PROJECT_ROOT *= "/.."   # Up to the effective root
include(PROJECT_ROOT * "/src/setup/simulations-setup.jl")
include(PROJECT_ROOT * "/src/modules/reciprocal-space-algorithm-general.jl")

using DelimitedFiles
using Dates

@doc raw"""
function RunHFRoutine(
    UU::Vector{Float64},
    VV::Vector{Float64},
    LL::Vector{Int64},
    δδ::Vector{Float64},
    ββ::Vector{Float64},
    p::Int64,
    Δm::Dict{String,Float64},
    Δn::Float64,
    g::Float64;
    FilePathOut::String=""
)

Returns: none if `FilePathOut` is specified, `ResultsTable::Matrix{Float64}` if
`FilePathOut` is unspecified.

`RunHFRoutine` takes as input `UU` (vector of local repulsions), `VV` (vector of
non local attractions), `LL` (vector of the square lattice dimensions), `δδ` 
(vector of dopings with respect to the half filled lattice), `ββ` (vector of 
inverse temperatures), `p` (maximum number of HF iterations), `Δm` (tolerance on
each order parameter) and `Δn` (tolerance on density in chemical potential
estimation), `g` (mixing parameter). It performs an iterative HF analysis over a
sequence of 2D square lattices for all the possible combinations of the 
specified parameters. Check the source files in the `/modules` folder for more
informations on the algorithm.
"""
function RunHFRoutine(
    Syms::Vector{String},       # Gap function symmetries
    UU::Vector{Float64},        # Local repulsion
    VV::Vector{Float64},		# Non-local attraction
    LL::Vector{Int64},          # Lattice size
    δδ::Vector{Float64},        # Doping
    ββ::Vector{Float64},        # Inverse temperature
    p::Int64,                   # Number of iterations
    Δm::Dict{String,Float64},   # Tolerance on magnetization
    Δn::Float64,                # Tolerance on density
    g::Float64;                 # Mixing parameter
    FilePathOut::String=""      # Output file
)
    
    if FilePathOut != ""
        Header = "# U, V, Lx, β, δ, Δ, Q, ΔT [calculated @ $(now())]\n"
        write(FilePathOut, Header)
    end

    SymStr = ""
    for Sym in Syms
        SymStr *= "$(Sym)-"
    end
    # Pop last "-" character
    SymStr = SymStr[1:end-1]

	Iterations = length(UU) * length(VV) * length(LL) * length(ββ) * length(δδ)
    i = 1    
    for (u,U) in enumerate(UU), 
        (v,V) in enumerate(VV), 
        Lx in LL, 
        β in ββ,
        δ in δδ

        L = [Lx, Lx]
        printstyled(
            "\e[2K\e[1GRun (i/$Iterations): $(SymStr)-wave HF at U=$U, V=$V, L=$Lx, β=$β, δ=$δ", 
            color=:yellow
        )
#        ResultsTable[i,1:5] .= [U,V,Lx,β,δ]
        ResultsVector::Matrix{Any} = [0 0]
        ResultsVector = hcat(ResultsVector, [U V Lx β δ])
        HFResults = RunHFAlgorithm(Syms,U,V,L,0.5+δ,β,p,Δm,Δn,g)
        
        m = Dict([Sym => HFResults[1][Sym] for Sym in Syms])
        Qs = Dict([Sym => HFResults[2][Sym] for Sym in Syms])        
		ResultsVector = hcat(ResultsVector[:,3:end], [m Qs HFResults[3]])

        i += 1

        if FilePathOut != ""
            open(FilePathOut, "a") do io
                writedlm(io, ResultsVector, ';')
        	end
        end
    end
    
    printstyled(
        "\e[2K\e[1GDone! Data saved at " * FilePathOut * "\n", color=:green
    )

end

# Main run
function main()
    DirPathOut = PROJECT_ROOT * "/simulations/Setup=$(Setup)/"
    mkpath(DirPathOut)

	HomeMadeRun = true
	if HomeMadeRun
		FilePathOut = DirPathOut * "$(HMSymsStr)-wave.txt"
		RunHFRoutine(HMSyms,UU,VV,LL,δδ,ββ,p,Δm,Δn,g;FilePathOut)
	elseif !HomeMadeRun
		AnisotropicSyms = ["px", "py", "d"]    
		for ASym in AnisotropicSyms
		    FilePathOut = DirPathOut * "$(ASym)-wave.txt"
		    RunHFRoutine([ASym],UU,VV,LL,δδ,ββ,p,Δm,Δn,g;FilePathOut)
		end

		FilePathOut = DirPathOut * "$Full-s-wave.txt"
		RunHFRoutine(["s", "s*"],UU,VV,LL,δδ,ββ,p,Δm,Δn,g;FilePathOut)
	end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
