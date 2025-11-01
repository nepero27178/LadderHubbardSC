#!/usr/bin/julia
using DelimitedFiles

@doc raw"""
function PlotOrderParameter(
	Phase::String,
    FilePathIn::String,
    DirPathOut::String;
    xVar::String=\"U\",
    pVar::String=\"V\",
    Skip::Int64=1,
    cs::Symbol=:imola50,
    RenormalizeHopping::Bool=true
)

Returns: none (plots saved at `DirPathOut`).

`PlotOrderParameter` takes as input `Phase` (string specifying the mean-field
phase, the allowed are \"AF\", \"SU/Singlet\", \"SU/Triplet\"), `FilePathIn`
(path to the data files), `DirPathOut` (path to the output directory). The
optional parameters are `xVar` and `pVar` (strings specifying respectively the x
variable and the parametric variable of the plot, the allowed are \"t\", \"U\",
\"V\", \"Lx\", \"δ\", \"β\"), `Skip` (integer indicating how many steps in the
parametric variable to be skipped between two plots), `cs` (colorscheme symbol).
The boolean option `RenormalizeHopping' allows for choosing to renormalize or 
not the hopping parameter.
"""
function PlotOrderParameter(
	Phase::String,						# Mean field phase
    FilePathIn::String,					# Data filepath
    DirPathOut::String;					# Output directory path
    xVar::String="U",                   # Specify x variable
    pVar::String="V",                   # Specify parametric variable
    Skip::Int64=1,                      # xVar skip parameter
    cs::Symbol=:imola50,                # Custom colorscheme
    RenormalizeHopping::Bool=true       # Conditional renormalization of t
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
	    "wp" => "w^{(\\bm{\\pi})}",
        # ...
    ])

    # Prepare title labels
    TitleLabels::Dict{String,String} = Dict([
	    "m" => "Magnetization",
	    "w0" => "\$w^{(\\mathbf{0})}\$",
	    "wp" => "\$w^{(\\bm{\\pi})}\$",
        # ...
    ])

    # Initialize directory structure
    DirPathOut *= "PlotOrderParameter/xVar=" * xVar * "/pVar=" * pVar * "/"
    mkpath(DirPathOut)
        
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
    uDF[xVar] = [NaN]                   # Remove from following cycle
    JJ::Vector = uDF[pVar][1:Skip:end]  # Parametric variable array
    uDF[pVar] = [NaN]                   # Remove from following cycle
    q = floor(Int64, length(colorschemes[cs]) / length(JJ) )

    # Check if the selected colorscheme is large enough
    if q==0
        @error "Your cs (ColorScheme) is not large enough. Select another."
        exit()
    end
    
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
        	TerminalMsg::String = "\e[2K\e[1GPlotting $(Phase) HF $(HF) data: "
        	FilePathOut::String = DirPathOut * "/" * HF
            rawTitle::String = TitleLabels[HF] * " ("

            AvoidList::Vector{String} = [xVar, pVar] #TODO Extend to pure hubbard
        	for Var in AllVars
        	    if !in(Var, AvoidList)
                    lVar = lDF[Var]
            	    TerminalMsg *= Var * "=$(lVar), "
            	    FilePathOut *= "_" * Var * "=$(lVar)"
            	    RS::String=""
            	    if !RenormalizeHopping && Var=="t"
            	        RS *= "\\tilde{t}="
            	    end
                    rawTitle *= "\$$(xVarLabels[Var])=" * RS * "$(lVar)\$, "
                    if Var=="Lx" # I am desperate about correct formatting
                        rawTitle = rawTitle[1:end-5] * "\$, "
                    elseif Var=="β" && β==Inf
                        rawTitle = rawTitle[1:end-6] * "\\infty\$, "
                    end
            	end
        	end
        	TerminalMsg = TerminalMsg[1:end-2] * " [x variable: " * xVar *
                ", parametric variable: " * pVar * "]"
        	FilePathOut *= ".pdf"
            rawTitle *= "varying \$" * xVarLabels[pVar] * "\$)"
        	printstyled("\e[2K\e[1G" * TerminalMsg, color=:yellow)
            
            # Cycle over parametric variable
            for (j,J) in enumerate(JJ)
                
                # Run local selection
            	lSelections = Selections .* (DF[pVar] .== J)
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
                    markercolor = colorschemes[cs][q*j],
                    markersize = 1.5,
                    linecolor = colorschemes[cs][q*j],
                    label = L"$%$(xVarLabels[pVar])=%$(J)$",
                    legendfonthalign = :left
                )
                title!(L"%$(rawTitle)")
            end
            
            # Save figure
            Plots.PGFPlotsX.push_preamble!(
                backend_object(P).the_plot, "\\usepackage{bm}"
            )
            savefig(P, FilePathOut)
		end
    end
    printstyled("\e[2K\e[1GDone! Plots saved at $(DirPathOut)\n", color=:green)
end

@doc raw"""
function PlotOrderParameter2D(
	Phase::String,
    FilePathIn::String,
    DirPathOut::String;
    xVar::String=\"U\",
    yVar::String=\"V\",
    cs::Symbol=:imola50
)

Returns: none (plots saved at `DirPathOut`).

`PlotOrderParameter2D` takes as input `Phase` (string specifying the mean-field
phase, the allowed are \"AF\", \"SU/Singlet\", \"SU/Triplet\"), `FilePathIn`
(path to the data files), `DirPathOut` (path to the output directory). The
optional parameters are `xVar` and `yVar` (strings specifying respectively the x
variable and the y variable of the plot, the allowed are \"t\", \"U\", \"V\", 
\"Lx\", \"δ\", \"β\"), `cs` (colorscheme symbol).
"""
function PlotOrderParameter2D(
	Phase::String,						# Mean field phase
    FilePathIn::String,					# Data filepath
    DirPathOut::String;					# Output directory path
    xVar::String="U",                   # Specify x variable
    yVar::String="V",                   # Specify y variable
    cs::Symbol=:imola50                 # Custom colorscheme
)
    
    FilePathOut = ""
    AllVars = ["t", "U", "V", "Lx", "δ", "β"]
    
    # Input safecheck
    if !in(xVar, AllVars)
        @error "Invalid x variable, choose one of $(AllVars)"
        exit()
    end
    
    if !in(yVar, AllVars)
        @error "Invalid parametric variable, choose one of $(AllVars)"
        exit()
    end
    
    if xVar==yVar
        @error "You have chosen xVar=yVar!"
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

    # Prepare title labels
    TitleLabels::Dict{String,String} = Dict([
	    "m" => "Magnetization",
	    "w0" => "\$w^{(\\mathbf{0})}\$",
	    "wp" => "\$w^{(\\bm{\\pi})}\$",
        # ...
    ])

    # Initialize directory structure
    DirPathOut *= "Heatmaps/xVar=" * xVar * "_yVar=" * yVar * "/"
    mkpath(DirPathOut)
        
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
    xx::Vector{Float64} = uDF[xVar]
    NumX::Int64 = length(xx)
    uDF[xVar] = [NaN]                   # Remove from following cycle
    yy::Vector{Float64} = uDF[yVar]
    NumY::Int64 = length(yy)
    uDF[yVar] = [NaN]                   # Remove from following cycle
    
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
        	    if !in(Var, [xVar, yVar])
        	        Selections = Selections .* (DF[Var] .== lDF[Var])
        	    end
        	end
    
            # Initialize plot    	
        	H = plot(
                size = (600,400),
                xlabel = L"$%$(xVarLabels[xVar])$",
                ylabel = L"$%$(xVarLabels[yVar])$",
                legend = :outertopright
            )
        	
        	# Write terminal message and file name
        	TerminalMsg::String = "\e[2K\e[1GPlotting $(Phase) HF $(HF) data: "
        	FilePathOut::String = DirPathOut * "/" * HF
            rawTitle::String = TitleLabels[HF] * " ("
        	for Var in AllVars
        	    if !in(Var, [xVar, yVar])
                    lVar = lDF[Var]
            	    TerminalMsg *= Var * "=$(lVar), "
            	    FilePathOut *= "_" * Var * "=$(lVar)"
                    rawTitle *= "\$$(xVarLabels[Var])=$(lVar)\$, "
                    if Var=="Lx" # I am desperate about correct formatting
                        rawTitle = rawTitle[1:end-5] * "\$, "
                    elseif Var=="β" && β==Inf
                        rawTitle = rawTitle[1:end-6] * "\\infty\$, "
                    end
            	end
        	end
        	TerminalMsg = TerminalMsg[1:end-2] * " [x variable: " * xVar *
                ", y variable: " * yVar * "]"
        	FilePathOut *= ".pdf"
            rawTitle = rawTitle[1:end-2] * ")"
        	printstyled("\e[2K\e[1G" * TerminalMsg, color=:yellow)
                             	
    	    # Define x, y variables and h
        	vv = DF["v"][Selections]
            hh::Matrix{Float64} = zeros(NumY,NumX)
            for j in 1:NumX
                hh[:,j] .= [eval(Meta.parse(
            	    vv[(j-1) * NumY + i]
                ))[HF] for i in 1:NumY]
            end

            # Plot parametrically
	        heatmap!(
                xx, yy, hh,
                color=cs
            )
            title!(L"%$(rawTitle)")

            # Save figure
            Plots.PGFPlotsX.push_preamble!(
                backend_object(H).the_plot, "\\usepackage{bm}"
            )
            savefig(H, FilePathOut)
        end
    end
    printstyled("\e[2K\e[1GDone! Plots saved at $(DirPathOut)\n", color=:green)
end

@doc raw"""
function PlotRMPs(
	Phase::String,
    FilePathIn::String,
    DirPathOut::String;
    xVar::String=\"U\",
    yVar::String=\"V\",
    cs::Symbol=:imola50
)

Returns: none (plots saved at `DirPathOut`).

`PlotRMPs` (Renormalized Model Parameters) takes as input `Phase` (string 
specifying the mean-field phase, the allowed are \"AF\", \"SU/Singlet\",
\"SU/Triplet\"), `FilePathIn` (path to the data files), `DirPathOut` (path to 
the output directory). The optional parameters are `xVar` and `yVar` (strings 
specifying respectively the x variable and the y variable of the plot, the 
allowed are \"t\", \"U\", \"V\", \"Lx\", \"δ\", \"β\"), `cs` (colorscheme 
symbol).
"""
function PlotRMPs(
	Phase::String,						# Mean field phase
    FilePathIn::String,					# Data filepath
    DirPathOut::String;					# Output directory path
    xVar::String="U",                   # Specify x variable
    yVar::String="V",                   # Specify y variable
    cs::Symbol=:imola50                 # Custom colorscheme
)
    
    FilePathOut = ""
    AllVars = ["t", "U", "V", "Lx", "δ", "β"]
    
    # Input safecheck
    if !in(xVar, AllVars)
        @error "Invalid x variable, choose one of $(AllVars)"
        exit()
    end
    
    if !in(yVar, AllVars)
        @error "Invalid parametric variable, choose one of $(AllVars)"
        exit()
    end
    
    if xVar==yVar
        @error "You have chosen xVar=yVar!"
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

    # Prepare title labels
    TitleLabels::Dict{String,String} = Dict([
	    "reΔ_tilde" => "\$\\mathrm{Re}" * 
            "\\lbrace\\tilde{\\Delta}_{\\mathbf{k}}\\rbrace\$",
	    "imΔ_tilde" => "\$\\mathrm{Im}" * 
            "\\lbrace\\tilde{\\Delta}_{\\mathbf{k}}\\rbrace\$",
	    "t_tilde" => "\$\\tilde{t}\$",
        # ...
    ])

    # Initialize directory structure
    DirPathOut *= "RMPs/xVar=" * xVar * "_yVar=" * yVar * "/"
    mkpath(DirPathOut)
        
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
    xx::Vector{Float64} = uDF[xVar]
    NumX::Int64 = length(xx)
    uDF[xVar] = [NaN]                   # Remove from following cycle
    yy::Vector{Float64} = uDF[yVar]
    NumY::Int64 = length(yy)
    uDF[yVar] = [NaN]                   # Remove from following cycle
    
    # List renormalized model parameters #TODO Extend to other phases
    ListRMPs::Vector{String} = [
        "reΔ_tilde", "imΔ_tilde", "t_tilde"
    ]

    # Cycle over renormalized model parameters
    for RMP in ListRMPs
    
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
        	    if !in(Var, [xVar, yVar])
        	        Selections = Selections .* (DF[Var] .== lDF[Var])
        	    end
        	end
    
            # Initialize plot    	
        	S = plot(
                size = (600,400),
                 xlabel = L"$%$(xVarLabels[xVar])$",
                 ylabel = L"$%$(xVarLabels[yVar])$",
                 zlabel = L"%$(TitleLabels[RMP])",
                legend = :outertopright
            )
        	
        	# Write terminal message and file name
        	TerminalMsg::String = "\e[2K\e[1GPlotting $(Phase) RMP $(RMP) " * 
                "data: "
        	FilePathOut::String = DirPathOut * "/" * RMP
            rawTitle::String = TitleLabels[RMP] * " ("
        	for Var in AllVars
        	    if !in(Var, [xVar, yVar])
                    lVar = lDF[Var]
            	    TerminalMsg *= Var * "=$(lVar), "
            	    FilePathOut *= "_" * Var * "=$(lVar)"
                    rawTitle *= "\$$(xVarLabels[Var])=$(lVar)\$, "
                    if Var=="Lx" # I am desperate about correct formatting
                        rawTitle = rawTitle[1:end-5] * "\$, "
                    elseif Var=="β" && β==Inf
                        rawTitle = rawTitle[1:end-6] * "\\infty\$, "
                    end
            	end
        	end
        	TerminalMsg = TerminalMsg[1:end-2] * " [x variable: " * xVar *
                ", y variable: " * yVar * "]"
        	FilePathOut *= ".pdf"
            rawTitle *= " \$k_\\ell=\\pi/3\$)"
        	printstyled("\e[2K\e[1G" * TerminalMsg, color=:yellow)
                             	
    	    # Define x, y variables and h
        	vv = DF["v"][Selections]
            hh::Matrix{Float64} = zeros(NumY,NumX)
            
            # Plot individually in order to adjust visual angles
            if RMP=="reΔ_tilde"
                
                for j in 1:NumX
                    hh[:,j] .= [eval(Meta.parse(
                	    vv[(j-1) * NumY + i]
                    ))["m"] for i in 1:NumY]
                end
                
                # Plot parametrically
                zz = hh .* (xx' .+ 8*yy)
                surface!(
                    xx, yy, zz,
                    color=cs,
                    label=L"%$(TitleLabels[RMP])",
                    camera=(30,25),
#                    zlim=(0.65,1),
#                    clim=(0.65,1)
                )
                title!(L"%$(rawTitle)")
            
            elseif RMP=="imΔ_tilde"
            
                for j in 1:NumX
                    hh[:,j] .= [eval(Meta.parse(
                	    vv[(j-1) * NumY + i]
                    ))["wp"] for i in 1:NumY]
                end
                
                # Plot parametrically
                zz = hh .* (2 .* yy * ones(1,length(xx)))
                surface!(
                    xx, yy, zz,
                    color=cs,
                    label=L"%$(TitleLabels[RMP])",
                    camera=(30,25),
#                    zlim=(0.65,1),
#                    clim=(0.65,1)
                )
                title!(L"%$(rawTitle)")
            
            elseif RMP=="t_tilde"
            
                for j in 1:NumX
                    hh[:,j] .= [eval(Meta.parse(
                	    vv[(j-1) * NumY + i]
                    ))["w0"] for i in 1:NumY]
                end

                # Plot parametrically
                zz = t*ones(size(hh)) .- hh .* (yy * ones(1,length(xx)))
	            surface!(
                    xx, yy, zz,
                    color=cs,
                    label=L"%$(TitleLabels[RMP])",
                    camera=(220,25),
                    zlim=(0.65,1),
                    clim=(0.65,1)
                )
                title!(L"%$(rawTitle)")
            end

            # Save figure
            Plots.PGFPlotsX.push_preamble!(
                backend_object(S).the_plot, "\\usepackage{bm}"
            )
            savefig(S, FilePathOut)
        end
    end
    printstyled("\e[2K\e[1GDone! Plots saved at $(DirPathOut)\n", color=:green)
end


@doc raw"""
function PlotRecord(
    Phase::String,
    DirPathIn::String,
    DirPathOut::String;
    rVar::String=\"g\",
    cs::Symbol=:imola25
)

Returns: none (plots saved at `DirPathOut`).

`PlotRecord` takes as input `Phase` (string specifying the mean-field phase, the
allowed are \"AF\", \"SU/Singlet\", \"SU/Triplet\"), `FilePathIn` (path to the
data files), `DirPathOut` (path to the output directory). The optional parameter 
is `rVar` (string specifying the recorded variable to plot), `cs` (colorscheme
symbol).
"""
function PlotRecord(
    Phase::String,                      # Mean field phase
    DirPathIn::String,                  # Data filepath
    DirPathOut::String;                 # Output directory path
    rVar::String="g",                   # Recorded variable
    cs::Symbol=:imola25                 # Custom colorscheme
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
        "wp" => "w^{(\\bm{\\pi})}",
        # ...
    ])

    # Prepare title labels
    TitleLabels::Dict{String,String} = Dict([
        "m" => "Magnetization",
        "w0" => "\$w^{(\\mathbf{0})}\$",
        "wp" => "\$w^{(\\bm{\\pi})}\$",
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
    	    "\$\\delta=$(δ)\$, "
    	    
        if β==Inf
    	    rawTitle *= "\$\\beta=\\infty\$)"
    	elseif β!=Inf
    	    rawTitle *= "\$\\beta=$(β)\$)"
    	end

        # Initialize plot
        P = plot(
            size = (600,400),
            xlim = (1,25),
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
            if length(yy)>100
                yy = yy[1:100]
                LineStyle = :dash
            end
            xx = [x for x in 1:length(yy)]
            plot!(
                xx,yy,
                markershape = :circle,
                markercolor = colorschemes[cs][j],
                markersize = 1.5,
                linestyle = LineStyle,
                linecolor = colorschemes[cs][j],
                label = L"$%$(pVarLabels[rVar])=%$(r)$",
            )

        end

        # Save figure
        Plots.PGFPlotsX.push_preamble!(
            backend_object(P).the_plot, "\\usepackage{bm}"
        )
        savefig(P, FilePathOut)
    end
    printstyled("\e[2K\e[1GDone! Plots saved at $(DirPathOut)\n", color=:green)
end
