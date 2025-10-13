#!/usr/bin/julia
using DelimitedFiles

@doc raw"""
...
"""
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

    # Prepare yVar labels
    yVarLabels::Dict{String,String} = Dict([
	    "m" => "m",
	    "w0" => "w^{(\\mathbf{0})}",
	    "wp" => "w^{(\\pi)}",
        # ...
    ])

    # Prepare title labels
    TitleLabels::Dict{String,String} = Dict([
	    "m" => "Magnetization",
	    "w0" => "\$w^{(\\mathbf{0})}\$",
	    "wp" => "\$w^{(\\pi)}\$",
        # ...
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
    
        # Cycle over simulated points
        for (w,t) in enumerate(uDF["t"]),
        	(u,U) in enumerate(uDF["U"]),
        	(v,V) in enumerate(uDF["V"]),
        	(l,L) in enumerate(uDF["Lx"]),
        	(d,δ) in enumerate(uDF["δ"]),
        	(b,β) in enumerate(uDF["β"])
        	
        	# Initialize local dataframe
        	lDF::Dict{String,Any} = Dict([
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
                ylabel = L"$%$(yVarLabels[HF])$",
                legend = :outertopright
            )
        	
        	# Write terminal message and file name
        	TerminalMsg::String = "\e[2K\e[1GPlotting $(Phase) HF $(HF) data\t"
        	FilePathOut::String = lDirPathOut * "/" * HF
            rawTitle::String = TitleLabels[HF] * " ("
        	for Var in AllVars
        	    if !in(Var, [xVar, pVar])
                    lVar = lDF[Var]
            	    TerminalMsg *= Var * "=$(lVar)\t"
            	    FilePathOut *= "_" * Var * "=$(lVar)"
                    rawTitle *= "\$$(xVarLabels[Var])=$(lVar)\$, "
                    if Var=="Lx" # I am desperate about correct formatting
                        rawTitle = rawTitle[1:end-5] * "\$, "
                    elseif Var=="β" && β==Inf
                        rawTitle = rawTitle[1:end-6] * "\\infty\$, "
                    end
            	end
        	end
        	TerminalMsg *= "x variable: " * xVar * 
                "\t\tparametric variable: " * pVar
        	FilePathOut *= ".pdf"
            rawTitle *= "varying " * xVarLabels[pVar] * ")"
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
                title!(L"%$(rawTitle)")
            end
            
            # Save figure
            savefig(P, FilePathOut)
		end
    end
    printstyled("\e[2K\e[1GDone! Plots saved at $(DirPathOut)\n", color=:green)
end

@doc raw"""
...
"""
function PlotRecord(
    Phase::String,
    DirPathIn::String,
    DirPathOut::String;
    rVar::String="g"                    # Recorded variable
)
        
    Files::Vector{String} = readdir(DirPathIn)
    rFiles::Dict{Float64,String} = Dict([])
    DataCols::Vector{String} = [""]

    # Check data structure and extract r values
    for (j,F) in enumerate(Files)
        FilePathIn::String = DirPathIn * "/" * F
        eqIndex::Int64 = findfirst('=', F)
        if F[1:eqIndex-1]==rVar

            # Read header
            Header = open(FilePathIn) do io
                readdlm(FilePathIn, ';', '\n')[1]
            end
            Start = findfirst('[',Header)
            Stop = findfirst(']',Header)
            if j==1
                DataCols = eval(Meta.parse( Header[Start:Stop] ))
            elseif j>1
                if DataCols != eval(Meta.parse( Header[Start:Stop] ))
                    @error "Bad data structure (DataCols inconsistent). " * 
                        "Check data integrity at " * DirPathIn
                    exit()
                end
            end
            r = parse(Float64, F[eqIndex+1:end-4]) # Remove .txt
            rFiles[r] = FilePathIn

        elseif F[1:eqIndex-1]!=rVar
            @error "Bad data structure (rVar not found). " * 
                "Check data integrity at " * DirPathIn
            exit()
        end
    end    

    # Prepare yVar labels
    yVarLabels::Dict{String,String} = Dict([
        "m" => "m",
        "w0" => "w^{(\\mathbf{0})}",
        "wp" => "w^{(\\pi)}",
        # ...
    ])

    # Prepare title labels
    TitleLabels::Dict{String,String} = Dict([
        "m" => "Magnetization",
        "w0" => "\$w^{(\\mathbf{0})}\$",
        "wp" => "\$w^{(\\pi)}\$",
        # ...
    ])

    pVarLabels::Dict{String,String} = Dict([
        "g" => "g",
        # ...
    ])
    
    for (hf,HF) in enumerate(DataCols)
        
        # Set output filepath
        FilePathOut::String = DirPathOut * "/" * HF * "_rVar=" * rVar * ".pdf"

        # Generate raw title
        rawTitle::String = TitleLabels[HF] * " (" *
            "\$t=$(t)\$, " *
            "\$U=$(U)\$, " *
            "\$V=$(V)\$, " *
            "\$L=$(L)\$, " *
    	    "\$\\delta=$(δ)\$, " *
    	    "\$\\beta=$(β)\$)"

        # Initialize plot
        P = plot(
            size = (600,400),
            xlabel = L"\mathrm{Step}",
            ylabel = L"$%$(yVarLabels[HF])$",
            legend = :outertopright
        )
        title!(L"%$(rawTitle)")

        # Cycle over data
        for (j,r) in enumerate( sort([key for key in keys(rFiles)]) )

            # Generate FilePathIn
            FilePathIn = rFiles[r]

            # Read and assign variables
            DataIn = open(FilePathIn) do io
                readdlm(FilePathIn, ';', comments=true, '\n')
            end

            yy = DataIn[:,hf]
            LineStyle = :solid
            if length(yy)>25
                yy = yy[1:25]
                LineStyle = :dash
            end
            xx = [x for x in 1:length(yy)]
            plot!(
                xx,yy,
                markershape = :circle,
                markercolor = ColorSchemes.imola25[j],
                markersize = 1.5,
                linestyle = LineStyle,
                linecolor = ColorSchemes.imola25[j],
                label = L"$%$(pVarLabels[rVar])=%$(r)$",
            )

        end

        # Save figure
        savefig(P, FilePathOut)
    end
    printstyled("\e[2K\e[1GDone! Plots saved at $(DirPathOut)\n", color=:green)
end
