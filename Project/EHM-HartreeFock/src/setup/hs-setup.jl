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

Phase = "AF" # Choose your phase
if !in(Phase, AllPhases)
	@error "Invalid phase, please modify at: " * SetupFilePath
	exit()
end

RenormalizeBands::Bool = false
if in(Phase,["AF","HybridSU-Singlet","HybridSU-Triplet"])
	RenormalizeBands = true
end


#-------------------------------------------------------------------------------
#------------------------------------ SYMS -------------------------------------
#-------------------------------------------------------------------------------

AllSingletSyms = ["s", "S", "d"]
AllTripletSyms = ["px", "py", "p+", "p-"]
Syms = ["s", "S"]
# Fake AF symmetry to respect general filestructure
occursin("AF",Phase) ? Syms=["π"] : 0

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

Setup = "Test[20]" # Choose your setup #TODO Use readline()
AvailableSetups = [
	"Test[20]",			# Test setup
	"A[250]",			# UV plane
	"B[250]",			# δV plane
	# OLD RUNS
	"A[128]",			# Pure Hubbard V=0
	"B[256]",			# Main run varying U, V, δ
	"-B[256]",			# Main run varying U, V, δ, but with -U
	"B[256]-t=0.7",		# Main run varying U, V, δ, but with rigid hopping renormalization
	"C[128]"	,			# Half filled model
	"D[256]",			# Sub-main run varying U, V, δ
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
elseif Setup=="Test[20]"
	# TEST-SETUP
	tt = [1.0]
	UU = [0.0]
	VV = [1.0, 2.0, 3.0]
	LL = [20]
	δδ = [0.0, 0.1, 0.2, 0.3]
	ββ = [100.0]
	p = 20
	Δv = TestΔv
	Δn = 1e-2
	g = 0.1


# --- MAIN A RUN ---
elseif Setup=="A[250]"
	tt = [1.0]
	UU = [U for U in 0.0:0.25:5.0]
	VV = [V for V in 0.0:0.25:5.0]
	LL = [150]
	δδ = [0.0, 0.2, 0.4]
	ββ = [100.0]
	p = 200
	Δv = MainΔv
	Δn = 1e-2
	g = 0.2

# --- MAIN B RUN ---
elseif Setup=="B[150]"
	tt = [1.0]
	UU = [0.0, 5.0, 10.0]
	VV = [V for V in 0.0:0.1:5.0]
	LL = [150]
	δδ = [δ for δ in 0.0:0.05:0.45]
	ββ = [100.0]
	p = 200
	Δv = MainΔv
	Δn = 1e-2
	g = 0.2

elseif Setup=="A[128]"
	tt = [2.0]
	UU = [U for U in 2.0:2.0:40.0]
	VV = [0.0]
	LL = [128]
	δδ = [δ for δ in 0.0:0.05:0.45]
	ββ = [100.0]
	p = 100
	Δv = MainΔv
	Δn = 1e-2
	g = 0.1
elseif Setup=="B[256]"
	tt = [1.0]
	UU = [4.0,12.0,20.0]
	VV = [V for V in 0.0:0.05:3.0]
	LL = [256]
	δδ = [δ for δ in 0.0:0.05:0.45]
	ββ = [100.0]
	p = 100
	Δv = MainΔv
	Δn = 1e-2
	g = 0.1
elseif Setup=="-B[256]"
	tt = [1.0]
	UU = [-4.0,-12.0,-20.0]
	VV = [V for V in 0.0:0.05:3.0]
	LL = [256]
	δδ = [δ for δ in 0.0:0.05:0.45]
	ββ = [100.0]
	p = 100
	Δv = MainΔv
	Δn = 1e-2
	g = 0.1
elseif Setup=="B[256]-t=0.7"
	tt = [0.7]
	UU = [4.0,12.0,20.0]
	VV = [V for V in 0.0:0.1:3.0]
	LL = [256]
	δδ = [δ for δ in 0.0:0.05:0.45]
	ββ = [100.0]
	p = 100
	Δv = MainΔv
	Δn = 1e-2
	g = 0.1
elseif Setup=="C[128]"
	tt = [1.0]
	UU = [U for U in 0.0:1.0:20.0]
	VV = [V for V in 0.0:0.5:10.0]
	LL = [128]
	δδ = [0.0]
	ββ = vcat([0.1, 1, 10], [β for β in 20:20:100])
	p = 100
	Δv = MainΔv
	Δn = 1e-2
	g = 0.1
elseif Setup=="D[256]"
	tt = [1.0]
	UU = [U for U in 0.0:0.5:15.0]
	VV = [V for V in 0.0:0.05:4.0]
	LL = [256]
	δδ = [0.0, 0.25, 0.45]
	ββ = [100.0]
	p = 1000
	Δv = MainΔv
	Δn = 1e-2
	g = 0.1
end
