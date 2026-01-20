#!/usr/bin/julia
# NOTE: This script is standalone importable and imports all simulations methods.

using LinearAlgebra
using Roots
using Random
LinearAlgebra.BLAS.set_num_threads(Threads.nthreads()) # Parallel optimization
using Integrals
using Elliptic
using DataFrames
using DelimitedFiles

PROJECT_METHODS_DIR = @__DIR__
include(PROJECT_METHODS_DIR * "/methods-physics.jl")
include(PROJECT_METHODS_DIR * "/methods-optimizations.jl")

@doc raw"""
function GetHFPs(
	Phase::String;
	Syms::Vector{String}=["s"]
)::Vector{String}

Returns: Hartree Fock Parameters labels for the given Phase.
"""
function GetHFPs(
	Phase::String;						# Mean field phase
	Syms::Vector{String}=["s"]			# Gap function symmetries
)::Vector{String}

	AF = false
	Singlet = false
	Triplet = false
	SymErr = "Invalid symmetries. $(Syms) is incoherent with $(Phase)."
	if in(Phase, ["AF", "FakeAF"])
		AF = true
	elseif in(Phase, ["SU-Singlet", "FakeSU-Singlet"])
		issubset(Syms, ["s", "S", "d"]) ? Singlet = true : throw(SymErr)
	elseif in(Phase, ["SU-Triplet", "FakeSU-Triplet"])
		issubset(Syms, ["px", "py", "p+", "p-"]) ? Triplet = true : throw(SymErr)
	end

	KeysList::Vector{String} = ["Empty"]
	AF ? KeysList = ["m", "w0", "wp"] : 0
	Singlet ? KeysList = vcat(["Δ$(Sym)" for Sym in Syms], "gS", "gd") : 0
	Triplet ? KeysList = vcat(["Δ$(Sym)" for Sym in Syms], "gS", "gd") : 0

	# Pure symmetry drop TODO Extension to Triplet
	if Singlet && (sort(Syms) == ["S", "s"] || Syms == ["d"])
		pop!(KeysList)
	end

	return KeysList

end

@doc raw"""
function FindRootμ(
	Phase::String,
	Parameters::Dict{String,Float64},
	K::Matrix{Vector{Float64}}
	v::Dict{String,Float64},
	nt::Float64,
	β::Float64;
	Δn::Float64=0.0,
	RenormalizeBands::Bool=true
)::Float64

Returns: a variational estimation for the chemical potential at density `nt`.

`FindRootμ` takes as input `Phase` (string specifying the mean-field phase, the
allowed are \"AF\", \"FakeAF\", \"SU-Singlet\", \"FakeSU-Singlet\",
\"SU-Triplet\", \"FakeSU-Triplet\"), `Parameters` (dictionary of model
parameters containing `t`, `U`, `V`), `K` (k-points in the BZ), `v`
(dictionary of real HF parameters), `nt` (target density) and `β` (inverse
temperature). It computes the root value for the chemical potential assuming
the matrix in reciprocal space to be the Nambu-Bogoliubov hamiltonian obtained
via the HF parameters `v`. The optimization is performed using the `Roots.jl`
library. The boolean option `RenormalizeBands' allows for choosing to
renormalize or not the free bands.
"""
function FindRootμ(
	Phase::String,						# Mean field phase
	Parameters::Dict{String,Float64},	# Model parameters t,U,V
	K::Matrix{Vector{Float64}},			# BZ grid
	v::Dict{String,Float64},			# HF parameters
	nt::Float64,						# Target density
	β::Float64;							# Inverse temperature
	Δn::Float64=0.0,					# Density tolerance
	debug::Bool=false,
	RenormalizeBands::Bool=true			# Conditional renormalization of t
)::Float64

	if nt < 0 || nt > 1
		@error "Invalid target density. Choose 0 ≤ nt ≤ 1." nt
		return
	end

	D::Int64 = 2 * prod(size(K))
	# Define function to be minimized
	δn(μ::Float64) = sum(
			GetKPopulation(Phase,Parameters,K,v,μ,β;RenormalizeBands)
		)/D - nt

	μ::Float64 = 0.0
	LowerBoundary::Float64 = 0.0
	UpperBoundary::Float64 = LowerBoundary
	if abs(δn(LowerBoundary)) > Δn
		if δn(LowerBoundary) > 0
			while δn(LowerBoundary) > 0
				if debug
					@warn "Moving down lower boundary"
				end
				LowerBoundary -= 1.0
			end
			UpperBoundary = LowerBoundary + 1.0
		elseif δn(UpperBoundary) < 0
			while δn(UpperBoundary) < 0
				if debug
					@warn "Moving up upper boundary"
				end
				UpperBoundary += 1.0
			end
			LowerBoundary = UpperBoundary - 1.0
		end
		μ = find_zero(δn, (LowerBoundary, UpperBoundary))
	end

	if debug
		n = sum( GetKPopulation(Phase,Parameters,K,v,μ,β) )/D
		@info "Optimal chemical potential and density:" μ n
	end

	return μ
end

@doc raw"""
function PerformHFStep(
	Phase::String,
	Parameters::Dict{String,Float64}
	K::Matrix{Vector{Float64}},
	v0::Dict{String,Float64},
	n::Float64,
	β::Float64;
	Syms::Vector{String}=[\"d\"],
	debug::Bool=false,
	RenormalizeBands::Bool=true
)::Tuple{Dict{String,Float64},Float64}

Returns: HF estimation for Cooper instability parameter based on input.

`PerformHFStep` takes as input `Phase` (string specifying the mean-field phase, 
the allowed are \"AF\", \"FakeAF\", \"SU-Singlet\", \"FakeSU-Singlet\",
\"SU-Triplet\", \"FakeSU-Triplet\"), `Parameters` (dictionary of model
parameters containing `t`, `U`, `V`), `K` (k-points in the BZ), `v0` (dictionary
of real HF initializers), `n` (density) and `β` (inverse temperature). It
performs a  single step of the iterative Hartree-Fock (HF) analysis starting
over a 2D square lattice whose dimensions are given by `L`. Among the optional
parameters is `Syms`, a strings vector specifying the superconducting gap
function symmetries to be simulated (allowed are: \"s\", \"S\", \"d\" for
the \"SU-Singlet\" phase, and \"px\", \"py\", \"p+\", \"p-\" for the
\"SU-Triplet\") phase. The function estimates `v` using as input the
Nambu-Bogoliubov hamiltonian obtained via the function `GetHamiltonian` and the
optimal chemical potential for the density `n` for this hamiltonian. The boolean
option `RenormalizeBands' allows for choosing to renormalize or not the
hopping parameter.
"""
function PerformHFStep(
	Phase::String,						# Mean field phase
	Parameters::Dict{String,Float64},	# Model parameters t,U,V
	K::Matrix{Vector{Float64}},			# BZ grid
	v0::Dict{String,Float64},			# HF initializers
	n::Float64,							# Density
	β::Float64;							# Inverse temperature
	Syms::Vector{String}=["s"],			# Gap function symmetries
	debug::Bool=false,					# Debug mode
	RenormalizeBands::Bool=true			# Conditional renormalization of t
)::Tuple{Dict{String,Float64},Float64}

	v = copy(v0)	
	LxLy = prod(size(K))	
	μ = FindRootμ(Phase,Parameters,K,v0,n,β;debug,RenormalizeBands)
	wk::Int64=0

	# Antiferromagnet
	if in(Phase, ["AF", "FakeAF"])
		m::Float64 = 0.0
		w0::Float64 = 0.0
		wpi::Float64 = 0.0
		
		t = Parameters["t"]
		if RenormalizeBands
			# Conditional renormalization of bands
			t -= v["w0"] * Parameters["V"]
		end
		
		for (i,q) in enumerate(K)
			wk = GetWeight(q; Sym="S-MBZ") # Avoid computational redundance
			k = q .* pi # Important: multiply k by pi
			if in(wk,[1,2,4]) # Allowed weights
				# Renormalized bands
				εk::Float64 = GetHoppingEnergy(t,k)

				# Renormalized gap
				reΔk::Float64 = v["m"] * (Parameters["U"] + 8*Parameters["V"])
				imΔk::Float64 = 2*v["wp"]*Parameters["V"] * StructureFactor("S",k)

				# Renormalized gapped bands
				Ek::Float64 = sqrt( εk^2 + reΔk^2 + imΔk^2 )

				# Fermi-Dirac factor
				FDk = FermiDirac(-Ek,μ,β) - FermiDirac(Ek,μ,β)

				if Ek!=0.0 # Otherwise add nothing
					m += wk * reΔk/Ek * FDk
					w0 -= wk * εk/Ek * FDk * StructureFactor("S",k)
					wpi += wk * imΔk/Ek * FDk * StructureFactor("S",k)
				end
			end
			wk = 0
		end
		
		v["m"] = m/(2*LxLy)
		v["w0"] = w0/(4*LxLy)
		v["wp"] = wpi/(4*LxLy)

	elseif in(Phase, ["SU-Singlet", "FakeSU-Singlet"])

		Δs::Float64 = 0.0 # s-wave
		ΔS::Float64 = 0.0 # s*-wave
		Δd::Float64 = 0.0 # d-wave
		gS::Float64 = 0.0 # s-wave part of g
		gd::Float64 = 0.0 # d-wave part of g

		t = Parameters["t"]
		if RenormalizeBands && "gS" in keys(v)
			# Conditional renormalization of bands
			t -= v["gS"]/2 * Parameters["V"]
		end

		for (i,q) in enumerate(K)
			wk = GetWeight(q) # Avoid computational redundance
			k = q .* pi # Important: multiply k by pi
			if in(wk,[1,2,4]) # Allowed weights
				# Free bands
				ξk = GetHoppingEnergy(t,k) - μ
				if RenormalizeBands && "gd" in keys(v)
					ξk += Parameters["V"] * v["gd"] * StructureFactor("d",k)
				end

				# Gap
				Δk::Float64 = 0.0
				for (key, value) in v
					if !in(key, ["gS", "gd"])
						key = String(key)
						key = String(chop(key, head=1, tail=0))
						Δk += value * StructureFactor(key,k)
					end
				end

				# Renormalized gapped bands
				Ek::Float64 = sqrt( ξk^2 + abs(Δk)^2 )

				# Gap factor
				sk::Float64 = 0.0
				ck::Float64 = 0.0
				if Ek!=0.0 # Otherwise add nothing
					sk = wk * Δk/Ek * tanh(β*Ek/2)
					ck = wk * ξk/Ek * tanh(β*Ek/2)
				end
				Δs -= sk
				ΔS += sk * StructureFactor("S",k)
				Δd += sk * StructureFactor("d",k)
				gS += (1-ck)/2 * StructureFactor("S",k)
				gd += (1-ck)/2 * StructureFactor("d",k)
			end
		end

		Δs *= Parameters["U"] / (2*LxLy)
		ΔS *= Parameters["V"] / LxLy
		Δd *= Parameters["V"] / LxLy
		gS /= LxLy
		gd /= LxLy

		"Δs" in keys(v0) ? v["Δs"] = Δs : 0
		"ΔS" in keys(v0) ? v["ΔS"] = ΔS : 0
		"Δd" in keys(v0) ? v["Δd"] = Δd : 0
		"gS" in keys(v0) ? v["gS"] = gS : 0
		"gd" in keys(v0) ? v["gd"] = gd : 0


	elseif Phase=="SU-Triplet"
		@error "Under construction"
		return
	end

	return v, μ
end

@doc raw"""
function RunHFAlgorithm(
	Phase::String,
	Parameters::Dict{String,Float64},
	L::Vector{Int64},
	n::Float64,
	β::Float64,
	p::Int64,
	Δv::Dict{String,Float64},
	Δn::Float64,
	g::Float64;
	v0::Dict{String,Float64}=Dict([]),
	Syms::Vector{String}=[\"d\"],
	verbose::Bool=false,
	debug::Bool=false,
	record::Bool=false,
	RenormalizeBands::Bool=true
)::Tuple{Dict{String,Float64}, Dict{String,Float64}, Float64, 
Dict{String,Vector{Float64}}}

Returns: Hartree-Fock (HF) estimation for `v`, sequential relative convergence 
parameters `Q` and computational time `ΔT`. If `record` is set to true, the
last output contains the entire evolution of the parameters.

`RunHFAlgorithm` takes as input `Phase` (string specifying the mean-field phase, 
the allowed are \"AF\", \"FakeAF\", \"SU-Singlet\", \"FakeSU-Singlet\",
\"SU-Triplet\", \"FakeSU-Triplet\"), `Parameters` (dictionary of model
parameters containing `t`, `U`, `V`), `L` (square lattice dimensions), `n`
(density), `β` (inverse temperature), `p` (maximum number of HF iterations),
`Δm` (tolerance on each order parameter) and `Δn` (tolerance on density in
chemical  potential estimation), `g` (mixing parameter). It performs an
iterative HF analysis over a 2D square lattice whose dimensions are given by
`L`. All of the function's positional arguments must be positive defined. The
algorithm is iterative and runs at most `p` times using at each reiteration the
result of the previous computation. The positional argument `v0` allows for
custom initialization of HF parameters. When left unspecified, `v0` elements
are randomly initialized. The `record` boolean positional parameter, when
activated, registers the entire evolution in the Vector{Float64} values of the
last output object. The boolean option `RenormalizeBands' allows for choosing
to renormalize or not the free bands.
"""
function RunHFAlgorithm(
	Phase::String,						# Mean field phase
	Parameters::Dict{String,Float64},	# Model parameters t,U,V
	L::Vector{Int64},					# [Lx, Ly]
	n::Float64,							# Density
	β::Float64,							# Inverse temperature
	p::Int64,							# Number of iterations
	Δv::Dict{String,Float64},			# Tolerance on each order parameter
	Δn::Float64,							# Tolerance on density difference
	g::Float64;							# Mixing parameter
	v0i::Dict=Dict([]),					# Initializers
	Syms::Vector{String}=["s"],			# Gap function symmetries
	verbose::Bool=false,
	debug::Bool=false,
	record::Bool=false,
	RenormalizeBands::Bool=true		# Conditional renormalization of t
)::Tuple{Dict{String,Any}, Dict{String,Any}}

	if verbose
		@info "Running HF convergence algorithm" Phase Syms Parameters L n β
		@info "Convergence parameters" p Δv Δn g
	end

	# Reciprocal space discretization (normalized to 1)
	Kx::Vector{Float64} = [kx for kx in -1:2/L[1]:1]
	popfirst!(Kx)
	Ky::Vector{Float64} = [ky for ky in -1:2/L[2]:1]
	popfirst!(Ky)
	K::Matrix{Vector{Float64}} = [ [kx,ky] for kx in Kx, ky in Ky ]

	# Get Hartree Fock Parameters labels
	HFPs = GetHFPs(Phase;Syms)

	# Initialize HF dictionaries
	v0::Dict{String,Float64} = Dict([])
	μ::Float64 = 0.0

	if v0i==Dict([])
		for HFP in HFPs
			v0[HFP] = rand()
		end
	elseif issubset(keys(v0i), HFPs)
		for HFP in HFPs
			v0[HFP] = copy(v0i[HFP])
		end
	end
	v = copy(v0) # Shallow copy of values
	Qs = copy(v0) # Copy NaN keys

	# Initialize record matrix
	Record::Dict{String,Vector{Float64}} = Dict([
		key => [ v0[key] ] for key in HFPs
	])

	# Recursive run
	i::Int64 = 1
	I::Int64 = p
	ΔT = @elapsed begin
		while i<=p

			if debug
				printstyled("\n---Step $i---\n", color=:yellow)
			end

			CurrentResults = PerformHFStep(
				Phase,
				Parameters,
				K,v0,n,β;
				Syms,
				debug,
				RenormalizeBands
			)
			v = copy(CurrentResults[1])
			μ = CurrentResults[2]
			for key in keys(v0)
				current = v[key]
				previous = v0[key]
				tolerance = Δv[key]
				Qs[key] = abs(current-previous) / tolerance
			end

			if all([Qs[key] for key in keys(v0)] .< 1)

				if verbose
					printstyled("\n---Converged at step $i---\n", color=:green)
				end
				I = i
				i = p+1

			elseif any([Qs[key] for key in keys(v0)] .>= 1)
				# g = rand() # Random mixing parameter
				for key in keys(v0)
					current = v[key]
					previous = v0[key]
					if record
						Record[key] = vcat(Record[key],g*current + (1-g)*previous)
					end
					v[key] = g*current + (1-g)*previous
				end

				if debug
					@info "Initializer and current step" v0 v
				end
				i += 1

			end
			v0 = copy(v)
		end
	end

	fMFT::Float64 = 0.0
	if all([Qs[key] for key in keys(v0)] .<= 1)

		if verbose
			@info "Algorithm has converged." v Qs
		end
		fMFT = FindRootμ(Phase,Parameters,K,v0,n,β;debug,RenormalizeBands)

	elseif any([Qs[key] for key in keys(v0)] .> 1)

		if verbose
			@info "Algorithm has not converged - v saved as NaN." v Qs Phase
		end

		# Substitute with NaN in order to plot blank points
		fMFT = NaN
		for key in keys(v0)
			v[key] = NaN
		end

	end

	Results::Dict{String,Any} = Dict([
		"HFPs" => v,
		"Record" => Record,
		"ChemicalPotential" => μ,
		"FreeEnergy" => fMFT
	])

	Performance::Dict{String,Any} = Dict([
		"Quality" => Qs,
		"Runtime" => ΔT,
		"Steps" => I
	])

	return Results, Performance
end
