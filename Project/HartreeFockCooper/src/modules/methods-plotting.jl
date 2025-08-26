#!/usr/bin/julia
using DelimitedFiles

function PlotVΔ(
    FilePathIn::String,
    DirPathOut::String
)

    DirPathOut *= "VΔ-Setup=$(Setup)/"
    mkpath(DirPathOut)

    DataIn = true
    open(FilePathIn) do io
        DataIn = readdlm(FilePathIn, ',', comments=true)
    end

    UU = unique(DataIn[:,1])
    VV = unique(DataIn[:,2])
    LL = unique(Int64.(DataIn[:,3]))
    ββ = unique(DataIn[:,4])
    δδ = unique(DataIn[:,5])

    # Cycler
    for (u,U) in enumerate(UU), (b,β) in enumerate(ββ), (l,L) in enumerate(LL)
        printstyled("\e[2K\e[1GPlotting HF d-wave data for β=$β", color=:yellow)
        FilePathOut = DirPathOut * "/U=$(U)_L=$(L)_β=$(β).pdf"
        P = plot(
            size = (600,400),
            xlabel = L"$V/t$",
            ylabel = L"$\Delta^{(d)}$",
            ylims = (0.0,0.5),
            legend = :topleft
        )
        if β==Inf
            title!(L"$d$-wave order parameter ($U=%$(U), L=%$(L), \beta=\infty$)")
        elseif β<Inf
            title!(L"$d$-wave order parameter ($U=%$(U), L=%$(L), \beta=%$(β)$)")
        end
        for (d,δ) in enumerate(δδ)
            
            Selections = (DataIn[:,1] .== U) .* (DataIn[:,3] .== L) .* 
            	(DataIn[:,4] .== β) .* (DataIn[:,5] .== δ) # Logical intersection
            VV = DataIn[Selections,2]
            (mm, QQ, ΔTΔT) = [DataIn[Selections,i] for i in 6:8]

            plot!(
                VV, mm,
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

function PlotδΔ(
    FilePathIn::String,
    DirPathOut::String
)

    DirPathOut *= "δΔ-Setup=$(Setup)/"
    mkpath(DirPathOut)

    DataIn = true
    open(FilePathIn) do io
        DataIn = readdlm(FilePathIn, ',', comments=true)
    end

    UU = unique(DataIn[:,1])
    VV = unique(DataIn[:,2])
    LL = unique(Int64.(DataIn[:,3]))
    ββ = unique(DataIn[:,4])
    δδ = unique(DataIn[:,5])

    # UU = UU[end-6:2:end]

    # Cycler
    for (u,U) in enumerate(UU), (b,β) in enumerate(ββ), (l,L) in enumerate(LL)
        printstyled("\e[2K\e[1GPlotting HF δΔ data for β=$β", color=:yellow)
        FilePathOut = DirPathOut * "/U=$(U)_L=$(L)_β=$(β).pdf"
        P = plot(
            size = (600,400),
            xlabel = L"$\delta \vphantom{V/t}$",
            ylabel = L"$\Delta^{(d)}$",
            ylims = (0.0,0.5),
            legend = :topright
        )
        if β==Inf
            title!(L"$d$-wave order parameter ($U=%$(U), L=%$(L), \beta=\infty$)")
        elseif β<Inf
            title!(L"$d$-wave order parameter ($U=%$(U), L=%$(L), \beta=%$(β)$)")
        end
        for (v,V) in enumerate(VV)
            
            Selections = (DataIn[:,2] .== V) .* (DataIn[:,3] .== L) .* 
                (DataIn[:,4] .== β) # Logical intersection
            δδ = DataIn[Selections,5]
            (mm, QQ, ΔTΔT) = [DataIn[Selections,i] for i in 6:8]

            plot!(
                δδ, mm,
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
