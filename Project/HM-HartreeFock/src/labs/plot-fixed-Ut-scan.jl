using DataFrames
using DelimitedFiles

# Includer
PROJECT_ROOT = @__DIR__
PROJECT_ROOT *= "/.."   # Up to the effective root
include(PROJECT_ROOT * "/setup/graphic-setup.jl")

function Plotm(
        DF::DataFrame
    )
    xx = unique(DF.x)
    δδ = unique(DF.δ)
    
    P = plot(
        size = (600,400),
        xlabel = L"x",
        ylabel = L"m",
        legend = :outertopright,
        legendfonthalign = :left
    )
    for (j,δ) in enumerate(δδ)
        df = DF[DF.δ.==δ,:]
        yy = [eval.(Meta.parse.(df.v[j]))["m"] for j in 1:length(df.x)]
        plot!(df.x, yy, color=ColorSchemes.lipari25[2*j], label=L"\delta=%$(δ)")
    end
    title!(L"Magnetization ($U=10$, $t=U/x$)")
    
    return P
end

function Plotμ(
        DF::DataFrame
    )
    xx = unique(DF.x)
    δδ = unique(DF.δ)
    
    P = plot(
        size = (600,400),
        xlabel = L"x",
        ylabel = L"\mu/t",
        legend = :outertopright,
        legendfonthalign = :left
    )
    for (j,δ) in enumerate(δδ)
        df = DF[DF.δ.==δ,:]
        yy = df.μ ./ df.t
        plot!(df.x, yy, color=ColorSchemes.lipari25[2*j], label=L"\delta=%$(δ)")
    end
    title!(L"Chemical potential ($U=10$, $t=U/x$)")
    
    return P
end

function PlotμΔ(
        DF::DataFrame
    )
    xx = unique(DF.x)
    δδ = unique(DF.δ)

    P = plot(
        size = (600,400),
        xlabel = L"x",
        ylabel = L"(\mu-\Delta)/t",
        legend = :outertopright,
        legendfonthalign = :left
    )
    for (j,δ) in enumerate(δδ)
        df = DF[DF.δ.==δ,:]
        ΔΔ = [eval.(Meta.parse.(df.v[j]))["m"] for j in 1:length(df.x)] .* 10
        yy = (df.μ .- ΔΔ) ./ df.t
        plot!(df.x, yy, color=ColorSchemes.lipari25[2*j], label=L"\delta=%$(δ)")
    end
    title!(L"Chemical potential eccedence ($U=10$, $t=U/x$)")

    return P
end

if abspath(PROGRAM_FILE) == @__FILE__
    FilePathIn = PROJECT_ROOT * "/../simulations/special/fixed-Ut-scan/vary-t/AF.txt"
    Data = readdlm(FilePathIn, ';')
    DF = DataFrame(Data[2:end,:], Data[1,:])
    P = Plotm(DF)
    savefig(P, "tmp.pdf")
end