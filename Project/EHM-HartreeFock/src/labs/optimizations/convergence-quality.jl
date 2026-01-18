using DataFrames
using DelimitedFiles

PROJECT_LABS_DIR = @__DIR__
include(PROJECT_LABS_DIR * "/../../modules/methods-simulating.jl")
include(PROJECT_LABS_DIR * "/../../modules/methods-IO.jl")

DF = ImportData(FilePathIn)
FilePathIn = "/home/nepero27178/Thesis/LadderHubbardSC/Project/EHM-HartreeFock/simulations/heatmap/Setup=Test[80]/Phase=FakeSU-Singlet/Syms=sS.txt"

