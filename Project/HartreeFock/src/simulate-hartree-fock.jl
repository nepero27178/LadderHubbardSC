#!/usr/bin/julia
PROJECT_ROOT = @__DIR__
PROJECT_ROOT *= "/.."   # Up to the effective root
include(PROJECT_ROOT * "/src/setup/simulations-setup.jl")
include(PROJECT_ROOT * "/src/modules/methods-simulating.jl")

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
    Phase::String,						# Mean field phase
    tt::Vector{Float64},                # Hopping amplitude
    UU::Vector{Float64},                # Local repulsion
    VV::Vector{Float64},	        	# Non-local attraction
    LL::Vector{Int64},                  # Lattice size
    δδ::Vector{Float64},                # Doping
    ββ::Vector{Float64},                # Inverse temperature
    p::Int64,                           # Number of iterations
    Δv::Dict{String,Float64},           # Tolerance on magnetization
    Δn::Float64,                        # Tolerance on density
    g::Float64;                         # Mixing parameter
    Syms::Vector{String}=["d"],		    # Gap function symmetries
    FilePathOut::String=""              # Output file
)

    # Phase discrimination
    keys = []
    if Phase=="AF"
		keys = ["m", "w0", "wp"]

	elseif Phase=="SU/Singlet"
		@error "Under construction"
		return

	elseif Phase=="SU/Triplet"
		@error "Under construction"
		return
	
	end
    
    # File initialization
    if FilePathOut != ""
        Header = "# [\"t\", \"U\", \"V\", \"Lx\", \"β\", \"δ\", " *
            "\"v\", \"Q\", \"ΔT\"] [calculated @ $(now())]\n"
        write(FilePathOut, Header)
    end

    # HF iterations
	Iterations = length(UU) * length(VV) * length(LL) * length(ββ) * length(δδ)
    i = 1    
    for t in tt,
        U in UU, 
        V in VV, 
        Lx in LL, 
        β in ββ,
        δ in δδ
        
        Parameters::Dict{String,Float64} = Dict([
            "t" => t,
            "U" => U,
            "V" => V
        ])

        L = [Lx, Lx]
        printstyled(
            "\e[2K\e[1GRun ($i/$Iterations): " *
			"$Phase HF at t=$t, U=$U, V=$V, L=$Lx, β=$β, δ=$δ", 
            color=:yellow
        )
        ResultsVector::Matrix{Any} = [0 0]  # Dummy placeholder
        ResultsVector = hcat(ResultsVector, [t U V Lx β δ])
        HFResults = RunHFAlgorithm(Phase,Parameters,L,0.5+δ,β,p,Δv,Δn,g)
        
        v = Dict([key => HFResults[1][key] for key in keys])
        Qs = Dict([key => HFResults[2][key] for key in keys])        
		ResultsVector = hcat(ResultsVector[:,3:end], [v Qs HFResults[3]])

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
    DirPathOut = PROJECT_ROOT * "/simulations/Phase=$(Phase)/Setup=$(Setup)/"
    mkpath(DirPathOut)
    FilePathOut = DirPathOut * Model * ".txt"
	RunHFRoutine(Phase,tt,UU,VV,LL,δδ,ββ,p,Δv,Δn,g;FilePathOut)
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
