#!/usr/bin/julia

# Arguments handler
if length(ARGS) != 1
    println("How to use this program?
Type the following: \$ julia ./plot.jl --mode
Where:
· mode = --scan / --heatmap / --RMPs / --record-g")
    exit()
else
    UserInput = ARGS
    Mode = UserInput[1][3:end]
end
InMode::String = Mode

# Includer
PROJECT_ROOT = @__DIR__
PROJECT_ROOT *= "/.."   # Up to the effective root
if in(Mode, ["scan", "heatmap", "record-g"])
    include(PROJECT_ROOT * "/src/setup/" * Mode * "-simulations-setup.jl")
elseif Mode=="RMPs" 
    include(PROJECT_ROOT * "/src/setup/heatmap-simulations-setup.jl")
    InMode = "heatmap"
else
    @error "Invalid argument. Use: mode = --scan / --heatmap / --RMPs " * 
        "/ --Record-g"
    exit()
end
include(PROJECT_ROOT * "/src/setup/graphic-setup.jl")
include(PROJECT_ROOT * "/src/modules/methods-plotting.jl")

function main()    
    DirPathIn = PROJECT_ROOT * "/simulations/Phase=" * Phase * "/" * 
        InMode * "/Setup=$(Setup)/"
    DirPathOut = PROJECT_ROOT * "/analysis/Phase=" * Phase * "/" * 
        Mode * "/Setup=$(Setup)/"
    mkpath(DirPathOut)
    if Mode=="scan"
        FilePathIn = DirPathIn * Model * ".txt"
        PlotOrderParameter(
            Phase,
            FilePathIn,
            DirPathOut;
            xVar="V",
            pVar="δ",
            cs=:winter,
            RenormalizeHopping
        )
        PlotOrderParameter(
            Phase,
            FilePathIn,
            DirPathOut;
            xVar="δ",          
            pVar="V",
            Skip=2,
            cs=:winter,
            RenormalizeHopping
        )
    elseif Mode=="heatmap"
        FilePathIn = DirPathIn * Model * ".txt"
        PlotOrderParameter2D(
            Phase,
            FilePathIn,
            DirPathOut;
            xVar="U",
            yVar="V",
            cs=:imola
        )
    elseif Mode=="RMPs"
        FilePathIn = DirPathIn * Model * ".txt"
        PlotRMPs(
            Phase,
            FilePathIn,
            DirPathOut;
            xVar="U",
            yVar="V",
            cs=:imola
        )
    elseif Mode=="record-g"
        PlotRecord(
            Phase,
            DirPathIn,
            DirPathOut;
            rVar="g"
        )
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
	main()
end
