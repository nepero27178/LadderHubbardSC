#!/usr/bin/julia
using DelimitedFiles
using Dates

# Arguments handler
if length(ARGS) != 1
    println("How to use this program?
Type the following: \$ julia ./simulate.jl --mode
Where:s
· --mode = --scan / --heatmap / --record-g")
    exit()
else
    UserInput = ARGS
    Mode = UserInput[1][3:end]
end

# Includer
PROJECT_ROOT = @__DIR__
PROJECT_ROOT *= "/.."   # Up to the effective root
if in(Mode, ["scan", "heatmap", "record-g"])
    include(PROJECT_ROOT * "/src/setup/" * Mode * "-simulations-setup.jl")
else
    @error "Invalid argument. Use: --mode = --scan / --heatmap / --record-g"
    exit()
end
include(PROJECT_ROOT * "/src/modules/methods-simulating.jl")

# Routines
@doc raw"""
function RunHFScan(
    UU::Vector{Float64},
    VV::Vector{Float64},
    LL::Vector{Int64},
    δδ::Vector{Float64},
    ββ::Vector{Float64},
    p::Int64,
    Δm::Dict{String,Float64},
    Δn::Float64,
    g::Float64;
    FilePathOut::String="",
    RenormalizeHopping::Bool=true
)

Returns: none if `FilePathOut` is specified, `ResultsTable::Matrix{Float64}` if
`FilePathOut` is unspecified.

`RunHFScan` takes as input `UU` (vector of local repulsions), `VV` (vector of
non local attractions), `LL` (vector of the square lattice dimensions), `δδ` 
(vector of dopings with respect to the half filled lattice), `ββ` (vector of 
inverse temperatures), `p` (maximum number of HF iterations), `Δm` (tolerance on
each order parameter) and `Δn` (tolerance on density in chemical potential
estimation), `g` (mixing parameter). It performs an iterative HF analysis over a
sequence of 2D square lattices for all the possible combinations of the 
specified parameters. Check the source files in the `/modules` folder for more
informations on the algorithm. The boolean option `RenormalizeHopping' allows 
for choosing to renormalize or not the hopping parameter.
"""
function RunHFScan(
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
    FilePathOut::String="",             # Output file
    RenormalizeHopping::Bool=true       # Conditional renormalization of t
)

    # Get Hartree Fock Parameters labels
    HFPs = GetHFPs(Phase)
    
    # File initialization
    if FilePathOut != ""
#=         Header = "# [\"t\", \"U\", \"V\", \"Lx\", \"β\", \"δ\", " *
            "\"v\", \"Q\", \"ΔT\"] [calculated @ $(now())]  \n" =#
            Header = "t;U;V;Lx;β;δ;v;Q;ΔT\n"
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
        
        # Run routine, all positional arguments here must be false
        HFResults = RunHFAlgorithm(
            Phase,Parameters,L,0.5+δ,β,
            p,Δv,Δn,g;
            RenormalizeHopping
        )
        
        v::Dict{String,Float64} = Dict([key => HFResults[1][key] for key in HFPs])
        Qs::Dict{String,Float64} = Dict([key => HFResults[2][key] for key in HFPs])        
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

@doc raw"""
...
"""
function RunHFRecord(
    Phase::String,						# Mean field phase
    t::Float64,                         # Hopping amplitude
    U::Float64,                         # Local repulsion
    V::Float64,                         # Non-local attraction
    Lx::Int64,                          # Lattice size
    δ::Float64,                         # Doping
    β::Float64,                         # Inverse temperature
    p::Int64,                           # Number of iterations
    Δv::Dict{String,Float64},           # Tolerance on magnetization
    Δn::Float64,                        # Tolerance on density
    gg::Vector{Float64};                # Mixing parameter
    Syms::Vector{String}=["d"],		    # Gap function symmetries
    DirPathOut::String="",              # Output file
    RenormalizeHopping::Bool=true       # Conditional renormalization of t
)::Dict{Float64,Dict{String,Vector{Float64}}}

    L = [Lx,Lx]

    # Get Hartree Fock Parameters labels
    HFPs = GetHFPs(Phase)
        
    # Initialize model parameters
    Parameters::Dict{String,Float64} = Dict([
        "t" => t,
        "U" => U,
        "V" => V
    ])

    printstyled(
		"Recording HF at t=$t\t U=$U\t V=$V\t L=$L\t β=$β\t δ=$δ\n", 
        color=:yellow
    )
    
    Record::Dict{Float64,Dict{String,Vector{Float64}}} = Dict([])
    for g in gg
        
        # Record routine
        HFResults = RunHFAlgorithm(
            Phase,Parameters,L,0.5+δ,β,
            p,Δv,Δn,g;
            verbose=true,
            record=true,
            RenormalizeHopping
        )
        ΔT::Float64 = HFResults[3]
        gRecord::Dict{String,Vector{Float64}} = Dict([
            key => HFResults[4][key] for key in HFPs
        ])

        # Write record on matrix
        RecordMatrix::Matrix{Float64} = zeros(
            length(values(
                gRecord[ HFPs[1] ]
            )),
            length( keys(HFPs) )
        )
        for (k,key) in enumerate(HFPs)
            RecordMatrix[:,k] = gRecord[key]
        end

        # Write on file
        if DirPathOut != ""
            FilePathOut = DirPathOut * "g=$(g).txt"

            # File initialization
            Header = "# $(HFPs) [calculated @ $(now())]\n"
            write(FilePathOut, Header)

            # Append recorded matrix
            open(FilePathOut, "a") do io
                writedlm(io, RecordMatrix, ';')
        	end
            printstyled(
                "\e[2K\e[1GDone! Data saved at " * FilePathOut *    
                "\n", color=:green
            )
        end
        Record[g] = gRecord
    end   
    return Record
end

# Main run
function main()
    # DirPathOut = PROJECT_ROOT * "/simulations/Phase=" * Phase * "/" * 
    DirPathOut = PROJECT_ROOT * "/simulations/" * 
        Mode * "/Setup=$(Setup)/"
    mkpath(DirPathOut)
    if in(Mode, ["scan", "heatmap"])
        # FilePathOut = DirPathOut * Model * ".txt"
	    FilePathOut = DirPathOut * Phase * ".txt"
	    RunHFScan(
	        Phase,
	        tt,UU,VV,
	        LL,δδ,ββ,
	        p,Δv,Δn,g;
	        FilePathOut,
	        RenormalizeHopping
        )
    elseif Mode=="record-g"
        RunHFRecord(
            Phase,
            t,U,V,
            L,δ,β,
            p,Δv,Δn,gg;
            DirPathOut,
            RenormalizeHopping
        )
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
