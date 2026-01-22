using GLMakie
using CairoMakie
using DataFrames
using DelimitedFiles

struct Simulation
	DF::DataFrame
	Setup::String
	Phase::String
	Syms::Vector{String}
end

struct GroupedPlot
	H::Figure
	DF::SubDataFrame
	FileName::String
end