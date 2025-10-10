#!/usr/bin/julia
using DelimitedFiles

function PlotVΔ(
	Syms::Vector{String},				# Gap function symmetries
    FilePathIn::String,					# Data filepath
    DirPathOut::String					# Output directory path
)

    DirPathOut *= "VΔ-Setup=$(Setup)/"
    mkpath(DirPathOut)

    DataIn = true
    open(FilePathIn) do io
        DataIn = readdlm(FilePathIn, ';', comments=true, '\n')
    end

    UU = unique(DataIn[:,1])
    VV = unique(DataIn[:,2])
    LL = unique(Int64.(DataIn[:,3]))
    ββ = unique(DataIn[:,4])
    δδ = unique(DataIn[:,5])

    # Cycler
    for Sym in Syms,
    	(u,U) in enumerate(UU),
    	(b,β) in enumerate(ββ),
    	(l,L) in enumerate(LL)

        printstyled(
        	"\e[2K\e[1GPlotting HF-δΔ $(Sym)-wave data for " *
			"U=$U, L=$L, β=$β", color=:yellow
        )
        FilePathOut = DirPathOut * "/$(Sym)-wave_U=$(U)_L=$(L)_β=$(β).pdf"
        P = plot(
            size = (600,400),
            xlabel = L"$V/t$",
            ylabel = L"$\Delta^{(%$(Sym))}$",
            # ylims = (0.0,0.25),
			xlims = (1.0,2.0),
			ylims = (0.0,0.04),
            legend = :topleft
        )
        if β==Inf
            title!(L"$%$(Sym)$-wave order parameter ($U/t=%$(U), L=%$(L), \beta=\infty$)")
        elseif β<Inf
            title!(L"$%$(Sym)$-wave order parameter ($U/t=%$(U), L=%$(L), \beta=%$(β)$)")
        end
        for (d,δ) in enumerate(δδ)
            
            Selections = (DataIn[:,1] .== U) .* # Logical intersection
            	(DataIn[:,3] .== L) .* 
            	(DataIn[:,4] .== β) .* 
            	(DataIn[:,5] .== δ)
            	
            VV = DataIn[Selections,2]
            (mm, QQ, ΔTΔT) = [DataIn[Selections,i] for i in 6:8]

			ΔΔ = [eval(
				Meta.parse(mm[i]) # Meta-programming: parse formatted string
			)[Sym] for i in 1:length(mm)]
            plot!(
                VV, ΔΔ,
                markershape = :circle,
                markercolor = TabColors[d],
                markersize = 1.5,
                linecolor = TabColors[d],
                label = L"$\delta=%$(δ)$",
                legendfonthalign = :left
            )
        end
        savefig(P, FilePathOut)
    end

    printstyled("\e[2K\e[1GDone! Plots saved at $(DirPathOut)\n", color=:green)

end

@doc raw"""
function PlotδΔ(
	Sym::Vector{String},
    FilePathIn::String,
    DirPathOut::String
)

Returns: none.

Under construction!
"""
function PlotδΔ(
	Syms::Vector{String},				# Gap function symmetries
    FilePathIn::String,					# Data filepath
    DirPathOut::String					# Output directory path
)

    DirPathOut *= "dΔ-Setup=$(Setup)/"
    mkpath(DirPathOut)

    DataIn = true
    open(FilePathIn) do io
        DataIn = readdlm(FilePathIn, ';', comments=true, '\n')
    end

    UU = unique(DataIn[:,1])
    VV = unique(DataIn[:,2])
    LL = unique(Int64.(DataIn[:,3]))
    ββ = unique(DataIn[:,4])
    δδ = unique(DataIn[:,5])

    # UU = UU[end-6:2:end]

    # Cycler
    for Sym in Syms,
    	(u,U) in enumerate(UU),
    	(b,β) in enumerate(ββ),
    	(l,L) in enumerate(LL)
    	
        printstyled(
        	"\e[2K\e[1GPlotting HF-δΔ $(Sym)-wave data for " *
			"U=$U, L=$L, β=$β", color=:yellow
        )
        FilePathOut = DirPathOut * "/$(Sym)-wave_U=$(U)_L=$(L)_β=$(β).pdf"
        P = plot(
            size = (600,400),
            xlabel = L"$\delta \vphantom{V/t}$",
            ylabel = L"$\Delta^{(%$(Sym))}$",
            # ylims = (0.0,0.25),
			xlims = (0.3,0.4),
			ylims = (0.0,0.06),
            legend = :topright
        )
        if β==Inf
            title!(L"$%$(Sym)$-wave order parameter ($U/t=%$(U), L=%$(L), \beta=\infty$)")
        elseif β<Inf
            title!(L"$%$(Sym)$-wave order parameter ($U/t=%$(U), L=%$(L), \beta=%$(β)$)")
        end
        for (v,V) in enumerate(VV)
            
            Selections = (DataIn[:,1] .== U) .*
				(DataIn[:,2] .== V) .*
				(DataIn[:,3] .== L) .* 
                (DataIn[:,4] .== β) # Logical intersection
            δδ = DataIn[Selections,5]
            (mm, QQ, ΔTΔT) = [DataIn[Selections,i] for i in 6:8]

			ΔΔ = [eval(
				Meta.parse(mm[i])	# Meta-programming: parse formatted string
			)[Sym] for i in 1:length(mm)]
            plot!(
                δδ, ΔΔ,
                markershape = :circle,
                markercolor = TabColors[v],
                markersize = 1.5,
                linecolor = TabColors[v],
                label = L"$V/t=%$(V)$",
                legendfonthalign = :left
            )
        end
        savefig(P, FilePathOut)
    end

    printstyled("\e[2K\e[1GDone! Plots saved at $(DirPathOut)\n", color=:green)

end

function PlotHeatmapVdΔ(
	Syms::Vector{String},				# Gap function symmetries
    FilePathIn::String,					# Data filepath
    DirPathOut::String					# Output directory path
)

	DirPathOut *= "Heatmap-VdΔ-Setup=$(Setup)/"
    mkpath(DirPathOut)

    DataIn = true
    open(FilePathIn) do io
        DataIn = readdlm(FilePathIn, ';', comments=true, '\n')
    end

    UU = unique(DataIn[:,1])
    VV = unique(DataIn[:,2])
    LL = unique(Int64.(DataIn[:,3]))
    ββ = unique(DataIn[:,4])
    δδ = unique(DataIn[:,5])

	# Cycler
    for Sym in Syms,
		(u,U) in enumerate([UU[1], UU[end]]),
    	(b,β) in enumerate([ββ[1], ββ[end]]),
    	(l,L) in enumerate([LL[1], LL[end]])
		
		FilePathOut = DirPathOut * "/$(Sym)-wave_U=$(U)_L=$(L)_β=$(β).pdf"
		Selections = (DataIn[:,1] .== U) .*
			(DataIn[:,3] .== L) .* 
            (DataIn[:,4] .== β) # Logical intersection

		NumV = length(VV)
    	Numδ = length(δδ)
		OrderParameters = zeros(Numδ,NumV)
		
		(mm, QQ, ΔTΔT) = [DataIn[Selections,i] for i in 6:8]
		ΔΔ = [eval(
			Meta.parse(mm[i])	# Meta-programming: parse formatted string
		)[Sym] for i in 1:length(mm)]

	    for vv in 1:NumV
	        OrderParameters[:,vv] = ΔΔ[ Numδ*(vv-1)+1 : Numδ*vv ]
		end

		heatmap(
			unique(VV), unique(δδ), OrderParameters, 
            xlabel=L"V",
            ylabel=L"$\delta$"
		)
		if β==Inf
            title!(L"$%$(Sym)$-wave order parameter ($U/t=%$(U), L=%$(L), \beta=\infty$)")
        elseif β<Inf
            title!(L"$%$(Sym)$-wave order parameter ($U/t=%$(U), L=%$(L), \beta=%$(β)$)")
        end
        savefig(FilePathOut)
	end
	printstyled("\e[2K\e[1GDone! Plots saved at $(DirPathOut)\n", color=:green)

end
