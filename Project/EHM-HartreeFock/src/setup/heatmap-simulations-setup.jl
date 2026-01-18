#!/usr/bin/julia
SetupFilePath = @__FILE__

#-------------------------------------------------------------------------------
#------------------------------------ PHASE ------------------------------------
#-------------------------------------------------------------------------------

AllPhases = [
	"AF",				# Renormalized antiFerromagnet
	"FakeAF",			# Renormalized antiFerromagnet, but with pure hopping
	"SU-Singlet",		# Singlet superconductor
	"FakeSU-Singlet",	# Singlet superconductor, but with pure hopping
	"SU-Triplet",		# Triplet superconductor
	"FakeSU-Triplet",	# Triplet superconductor, but with pure hopping
	# "HybridSU"			# Hybrid (singlet+triplet) superconductor
]

Phase = "SU-Singlet" # Choose your phase
if !in(Phase, AllPhases)
	@error "Invalid phase, please modify at: " * SetupFilePath
	exit()
end

RenormalizeBands::Bool = false #TODO => RenormalizeBands::Bool
if in(Phase,["AF","SU-Singlet","SU-Triplet"])
	RenormalizeBands = true
end

#-------------------------------------------------------------------------------
#------------------------------------ SYMS -------------------------------------
#-------------------------------------------------------------------------------

AllSingletSyms = ["s", "S", "d"]
AllTripletSyms = ["px", "py", "p+", "p-"]
Syms = ["s", "S"]
if in(Phase, ["SU-Singlet", "FakeSU-Singlet"])
	if !issubset(Syms, AllSingletSyms)
		@error "Invalid symmetries. $(Syms) is incoherent with $(Phase)." *
			"Please modify at: " * SetupFilePath
	end
elseif in(Phase, ["SU-Triplet", "FakeSU-Triplet"])
	if !issubset(Syms, AllTripletSyms)
		@error "Invalid symmetries. $(Syms) is incoherent with $(Phase)." *
			"Please modify at: " * SetupFilePath
	end
end

#-------------------------------------------------------------------------------
#------------------------------------ SETUP ------------------------------------
#-------------------------------------------------------------------------------

Setup = "A[128]" # Choose your setup #TODO Use readline()
AvailableSetups = [
	"Test[80]",			# Test setup
	"A[128]",			# UV plane
	"-A[128]",			# UV plane, but with -U
	"B[128]",			# δV plane
	"-B[128]",			# δV plane, but with -U
	"B[128]-t=0.7",		# δV plane, but with rigid hopping renormalization
	"C[128]",			# δV plane
	"-C[128]",			# δV plane, but with -U
	"D[128]",
	"D[128]-Zoom"
]

TestΔv::Dict{String,Float64} = Dict([
	key => 1e-3 for key in [
		"m","w0","wp",
		"Δs","ΔS","Δd",
		"gS","gd",
		"Δpx","Δpy","Δp+","Δp-"
	]
])
MainΔv::Dict{String,Float64} = Dict([
	key => 1e-4 for key in [
		"m","w0","wp",
		"Δs","ΔS","Δd",
		"gS","gd",
		"Δpx","Δpy","Δp+","Δp-"
	]
])

if !in(Setup, AvailableSetups)
	@error "Invalid setup, please modify at: " * SetupFilePath
elseif Setup=="Test[80]"
	# TEST-SETUP
	tt = [1.0]
	UU = [4.0]
	VV = [V for V in 0.0:0.5:4.0]
	LL = [80]
	δδ = [δ for δ in 0.0:0.05:0.45]
	ββ = [100.0]
	p = 100
	Δv = TestΔv
	Δn = 1e-2
	g = 0.5
elseif Setup=="A[128]"
	tt = [1.0]
	UU = [U for U in 0.0:1.0:20.0]
	VV = [V for V in 0.0:0.1:4.0]
	LL = [128]
	δδ = [0.2]
	ββ = [100.0]
	p = 200
	Δv = MainΔv
	Δn = 1e-2
	g = 0.5
elseif Setup=="-A[128]"
	tt = [1.0]
	UU = [U for U in -20.0:1.0:10.0]
	VV = [V for V in 0.0:0.1:4.0]
	LL = [128]
	δδ = [0.2]
	ββ = [100.0]
	p = 200
	Δv = MainΔv
	Δn = 1e-2
	g = 0.5
elseif Setup=="B[128]"
	tt = [1.0]
	UU = [4.0]
	VV = [V for V in 0.0:0.1:4.0]
	LL = [128]
	δδ = [δ for δ in 0.0:0.01:0.49]
	ββ = [100.0]
	p = 200
	Δv = MainΔv
	Δn = 1e-2
	g = 0.5
elseif Setup=="-B[128]"
	tt = [1.0]
	UU = [-4.0]
	VV = [V for V in 0.0:0.1:4.0]
	LL = [128]
	δδ = [δ for δ in 0.0:0.01:0.49]
	ββ = [100.0]
	p = 200
	Δv = MainΔv
	Δn = 1e-2
	g = 0.5
elseif Setup=="B[128]-t=0.7"
	tt = [0.7]
	UU = [4.0]
	VV = [V for V in 0.0:0.1:4.0]
	LL = [128]
	δδ = [δ for δ in 0.0:0.01:0.49]
	ββ = [100.0]
	p = 200
	Δv = MainΔv
	Δn = 1e-2
	g = 0.5
elseif Setup=="C[128]"
	tt = [1.0]
	UU = [12.0]
	VV = [V for V in 0.0:0.1:4.0]
	LL = [128]
	δδ = [δ for δ in 0.0:0.01:0.49]
	ββ = [100.0]
	p = 200
	Δv = MainΔv
	Δn = 1e-2
	g = 0.5
elseif Setup=="-C[128]"
	tt = [1.0]
	UU = [-12.0]
	VV = [V for V in 0.0:0.1:4.0]
	LL = [128]
	δδ = [δ for δ in 0.0:0.01:0.49]
	ββ = [100.0]
	p = 200
	Δv = MainΔv
	Δn = 1e-2
	g = 0.5
elseif Setup=="D[128]"
	tt = [1.0]
	UU = [U for U in 0.0:1.0:20.0]
	VV = [V for V in 0.0:0.5:10.0]
	LL = [128]
	δδ = [0.0]
	ββ = [100.0, 50.0, 10.0]
	p = 200
	Δv = MainΔv
	Δn = 1e-2
	g = 0.5
elseif Setup=="D[128]-Zoom"
	tt = [1.0]
	UU = [U for U in 0.0:0.25:15]
	VV = [V for V in 0.0:0.1:2.5]
	LL = [128]
	δδ = [0.0]
	ββ = [100.0, 10.0]
	p = 200
	Δv = MainΔv
	Δn = 1e-2
	g = 0.5
end
