#!/usr/bin/julia

# Arguments handler
if length(ARGS) != 1
    println("How to use this program?
Type the following: \$ julia ./ising2D_metro.jl mode
Where:
· mode = \"Scan\" / \"Record_g\"")
    exit()
else
    UserInput = ARGS
    Mode = UserInput[1]
end

# Includer
PROJECT_ROOT = @__DIR__
PROJECT_ROOT *= "/.."   # Up to the effective root
if Mode=="Scan"
    include(PROJECT_ROOT * "/src/setup/scan-simulations-setup.jl")
elseif Mode=="Record_g"
    include(PROJECT_ROOT * "/src/setup/record-g-simulations-setup.jl")
else
    @error "Invalid argument. Use: mode = Scan / Record_g"
    exit()
end
include(PROJECT_ROOT * "/src/setup/graphic-setup.jl")
include(PROJECT_ROOT * "/src/modules/methods-plotting.jl")

function main()    
    DirPathIn = PROJECT_ROOT * "/simulations/Phase=" * Phase * "/" * 
        Mode * "/Setup=$(Setup)/"
    DirPathOut = PROJECT_ROOT * "/analysis/Phase=" * Phase * "/" * 
        Mode * "/Setup=$(Setup)/"
    mkpath(DirPathOut)
    if Mode=="Scan"
        FilePathIn = DirPathIn * Model * ".txt"
        PlotOrderParameter(
            Phase,
            FilePathIn,
            DirPathOut;
            xVar="V",
            pVar="δ"
        )
        PlotOrderParameter(
            Phase,
            FilePathIn,
            DirPathOut;
            xVar="δ",          
            pVar="V"
        )
    elseif Mode=="Record_g"
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
