#!/usr/bin/julia
using DelimitedFiles

function PlotUm(
    FilePathIn::String,
    DirPathOut::String
)

    DirPathOut *= "uM-Setup=$(Setup)/"
    mkpath(DirPathOut)

    DataIn = true
    open(FilePathIn) do io
        DataIn = readdlm(FilePathIn, ',', comments=true)
    end

    UU = unique(DataIn[:,1])
    LL = unique(Int64.(DataIn[:,3]))
    ββ = unique(DataIn[:,4])
    δδ = unique(DataIn[:,5])

    # Cycler
    for (b,β) in enumerate(ββ), (l,L) in enumerate(LL)
        printstyled("\e[2K\e[1GPlotting HF mU data for β=$β", color=:yellow)
        FilePathOut = DirPathOut * "/L=$(L)_β=$(β).pdf"
        P = plot(
            size = (600,400),
            xlabel = L"$U/t$",
            ylabel = L"$m$",
            ylims = (0.0,0.5),
            legend = :topleft
        )
        if β==Inf
            title!(L"Magnetization ($L=%$(L), \beta=\infty$)")
        elseif β<Inf
            title!(L"Magnetization ($L=%$(L), \beta=%$(β)$)")
        end
        for (d,δ) in enumerate(δδ)
            
            Selections = (DataIn[:,4] .== β) .* (DataIn[:,5] .== δ)
            UU = DataIn[Selections,1]
            (mm, QQ, ΔTΔT) = [DataIn[Selections,i] for i in 6:8]

            plot!(
                UU, mm,
                markershape = :circle,
                markercolor = TabColors[d],
                markersize = 1.5,
                label = L"$\delta=%$(δ)$"
            )
        end
        savefig(P, FilePathOut)
    end

    printstyled("\e[2K\e[1GDone! Plots saved at $(DirPathOut)\n", color=:green)

end

function Plotδm(
    FilePathIn::String,
    DirPathOut::String
)

    DirPathOut *= "δM-Setup=$(Setup)/"
    mkpath(DirPathOut)

    DataIn = true
    open(FilePathIn) do io
        DataIn = readdlm(FilePathIn, ',', comments=true)
    end

    UU = unique(DataIn[:,1])
    LL = unique(Int64.(DataIn[:,3]))
    ββ = unique(DataIn[:,4])
    δδ = unique(DataIn[:,5])

    UU = UU[end-6:2:end]

    # Cycler
    for (b,β) in enumerate(ββ), (l,L) in enumerate(LL)
        printstyled("\e[2K\e[1GPlotting HF δM data for β=$β", color=:yellow)
        FilePathOut = DirPathOut * "/L=$(L)_β=$(β).pdf"
        P = plot(
            size = (600,400),
            xlabel = L"$\delta$",
            ylabel = L"$m$",
            ylims = (0.0,0.5),
            legend = :bottomleft
        )
        if β==Inf
            title!(L"Magnetization ($L=%$(L), \beta=\infty$)")
        elseif β<Inf
            title!(L"Magnetization ($L=%$(L), \beta=%$(β)$)")
        end
        for (u,U) in enumerate(UU)
            
            Selections = (DataIn[:,1] .== U) .* (DataIn[:,4] .== β)
            δδ = DataIn[Selections,5]
            (mm, QQ, ΔTΔT) = [DataIn[Selections,i] for i in 6:8]

            plot!(
                δδ, mm,
                markershape = :circle,
                markercolor = TabColors[u],
                markersize = 1.5,
                label = L"$U/t=%$(Int64(U))$"
            )
        end
        savefig(P, FilePathOut)
    end

    printstyled("\e[2K\e[1GDone! Plots saved at $(DirPathOut)\n", color=:green)

end
