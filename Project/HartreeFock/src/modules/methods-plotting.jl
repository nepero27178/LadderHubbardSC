#!/usr/bin/julia
using DelimitedFiles

function PlotOrderParameter(
	Phase::String,						# Mean field phase
    FilePathIn::String,					# Data filepath
    DirPathOut::String;					# Output directory path
    xVar::String="U",                   # Specify x variable
    pVar::String="V",                   # Specify parametric variable
)
    
    FilePathOut = ""
    AllVars = ["t", "U", "V", "Lx", "δ", "β"]
    
    # Input safecheck
    if !in(xVar, AllVars)
        @error "Invalid x variable, choose one of $(AllVars)"
        exit()
    end
    
    if !in(pVar, AllVars)
        @error "Invalid parametric variable, choose one of $(AllVars)"
        exit()
    end
    
    if xVar==pVar
        @error "You have chosen xVar=pVar!"
        exit()
    end

    # Prepare xVar labels
    xVarLabels::Dict{String,String} = Dict([
	"t" => "t",
	"U" => "U",
	"V" => "V",
	"Lx" => "L_x",
	"δ" => "\\delta",
	"β" => "\\beta"
    ])

    # Initialize directory structure
    DirPathOut *= "PlotOrderParameter/xVar=" * xVar * "/pVar=" * pVar * "/"
        
    # Read header
    Header = open(FilePathIn) do io
        readdlm(FilePathIn, ';', '\n')[1]
    end
    Start = findfirst('[',Header)
    Stop = findfirst(']',Header)
    DataCols = eval(Meta.parse( Header[Start:Stop] ))
    
    # Read and assign variables
    DataIn = open(FilePathIn) do io
        readdlm(FilePathIn, ';', comments=true, '\n')
    end
    
    # Create DFs
    DF::Dict{String,Any} = Dict([])
    uDF::Dict{String,Any} = Dict([])
    for (v,Var) in enumerate(DataCols)
        Data = DataIn[:,v]
        if Var=="Lx"
            Data=Int64.(Data)
        end
        DF[Var] = Data
        uDF[Var] = unique(Data)
    end
    uDF[xVar] = [NaN] # Remove from following cycle
    pp::Vector = uDF[pVar]
    uDF[pVar] = [NaN] # Remove from following cycle
    
    # List HF parameters
    ListHF::Vector{String} = [key for key in keys(
        eval(Meta.parse(
            DF["v"][1]
        ))
    )]
    
    # Cycle over HF parameters
    for HF in ListHF
    
        # Initialize local data structure
        lDirPathOut = DirPathOut * HF * "/"
        mkpath(lDirPathOut)
        LabelHF::String = HF
        if HF=="w0"
            LabelHF = "w^{(\\mathbf{0})}"
        elseif HF=="wp"
            LabelHF = "w^{(\\pi)}"
        end
    
        # Cycle over simulated points
        for (w,t) in enumerate(uDF["t"]),
        	(u,U) in enumerate(uDF["U"]),
        	(v,V) in enumerate(uDF["V"]),
        	(l,L) in enumerate(uDF["Lx"]),
        	(d,δ) in enumerate(uDF["δ"]),
        	(b,β) in enumerate(uDF["β"])
        	
        	# Initialize local dataframe
        	lDF::Dict{String,Float64} = Dict([
        	    "t" => t,
        	    "U" => U,
        	    "V" => V,
        	    "Lx" => L,
        	    "δ" => δ,
        	    "β" => β
        	])
        	
        	# Select entries
        	Selections::Vector{Bool} = [true for _ in 1:length(DF["t"])]
        	for Var in AllVars
        	    if !in(Var, [xVar, pVar])
        	        Selections = Selections .* (DF[Var] .== lDF[Var])
        	    end
        	end
    
            # Initialize plot    	
        	P = plot(
                size = (600,400),
                xlabel = L"$%$(xVarLabels[xVar])$",
                ylabel = L"$%$(LabelHF)$",
                legend = :outertopright
            )
        	
        	# Write terminal message and file name
        	TerminalMsg::String = "\e[2K\e[1GPlotting $(Phase) HF $(HF) data for "
        	FilePathOut::String = lDirPathOut * "/" * HF
        	for Var in AllVars
        	    if !in(Var, [xVar, pVar])
            	    TerminalMsg *= Var * "=$(lDF[Var]), "
            	    FilePathOut *= "_" * Var * "=$(lDF[Var])"
            	end
        	end
        	TerminalMsg *= "x variable: " * xVar * ", parametric variable: " * pVar
        	FilePathOut *= ".pdf"
        	printstyled("\e[2K\e[1G" * TerminalMsg, color=:yellow)
            
            # Cycle over parametric variable
            for (j,P) in enumerate(pp)
                
                # Run local selection
            	lSelections = Selections .* (DF[pVar] .== P)
            	vv = DF["v"][lSelections]
        	
        	    # Define x and y variables
            	xx::Vector{Float64} = DF[xVar][lSelections]
            	yy::Vector = [eval(Meta.parse(
            	    vv[i]
		        ))[HF] for i in 1:length(vv)]

                # Plot parametrically
		        plot!(
                    xx, yy,
                    markershape = :circle,
                    markercolor = ColorSchemes.imola25[j], # TabColors[j],
                    markersize = 1.5,
                    linecolor = ColorSchemes.imola25[j], # TabColors[j],
                    label = pVar * "=$(P)",
                    legendfonthalign = :left
                )
            end
            
            # Save figure
            savefig(P, FilePathOut)
		end
    end
    printstyled("\e[2K\e[1GDone! Plots saved at $(DirPathOut)\n", color=:green)
end
