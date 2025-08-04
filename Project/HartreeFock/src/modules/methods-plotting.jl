#!/usr/bin/julia
using DelimitedFiles

function PlotHFData(
    FilePathIn::String,
    DirPathOut::String;
    CustomLL::Vector{Any}=[],
    Customββ::Vector{Any}=[],
)

    DataIn = true
    open(FilePathIn) do io
        DataIn = readdlm(FilePathIn, ',', comments=true)
    end

    if length(CustomLL)>0
        LL = CustomLL
    elseif length(CustomLL)==0
        LL = unique(Int64.(DataIn[:,1]))
    end

    if length(Customββ)>0
        ββ = Customββ
    elseif length(Customββ)==0
        ββ = unique(DataIn[:,2])
    end

    for (b,β) in enumerate(ββ)
        printstyled("\e[2K\e[1GPlotting HF data for β=$β", color=:yellow)
        FilePathOut = DirPathOut * "/t=$(t)_β=$(β).pdf"
        P = plot(
            xlabel = L"$\delta = n-0.5$",
            ylabel = L"$m$",
            ylims = (-0.5,0.5)
        )
        for (l,L) in enumerate(LL)
            title!(L"Magnetization ($t/U=%$(t), \beta=%$(β)$)")
            Selections = (DataIn[:,1] .== L) .* (DataIn[:,2] .== β)
            (δδ, mm, QQ, ΔTΔT) = [DataIn[Selections,i] for i in 3:6]

            plot!(
                δδ, mm,
                markershape = :circle,
                markercolor = TabColors[l],
                label = L"$L=%$(L)$"
            )
        end
        savefig(P, FilePathOut)
    end

    printstyled("\e[2K\e[1GDone! Plots saved at $(DirPathOut)\n", color=:green)

end
