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
include(PROJECT_METHODS_DIR * "/methods-IO.jl")

@doc raw"""
function FindRootμ(
	Phase::String,
	Parameters::Dict{String,Float64},
	K::Matrix{Vector{Float64}}
	v::Dict{String,Float64},
	nt::Float64,
	β::Float64;
	Δn::Float64=0.0,
	μ0::Float64=0.0,
	RenormalizeBands::Bool=true,
	OptimizeBZ::Bool=true,
	debug::Bool=false
)::Float64

Returns: optimal chemical potential.
"""
function FindRootμ(
	Phase::String,						# Mean field phase
	Parameters::Dict{String,Float64},	# Model parameters t,U,V
	K::Matrix{Vector{Float64}},			# BZ grid
	v::Dict{String,Float64},				# HF parameters
	nt::Float64,							# Target density
	β::Float64;							# Inverse temperature
	Δn::Float64=1e-3,					# Density tolerance
	μ0::Float64=0.0,						# Initial guess
	RenormalizeBands::Bool=true,			# Conditional renormalization of t
	OptimizeBZ::Bool=true,				# Conditional BZ optimization
	debug::Bool=false,					# Debug mode
)::Float64

	if nt < 0 || nt > 1
		@error "Invalid target density. Choose 0 ≤ nt ≤ 1." nt
		return
	end

	μ::Float64 = 0.0
	if nt != 0.5 # Speed up at half-filling

		D::Int64 = 2 * prod(size(K))
		# Define function to be minimized
		δn(μ::Float64) = sum(
				GetKPopulation(Phase,Parameters,K,v,μ,β;RenormalizeBands,OptimizeBZ)
			)/D - nt

		LowerBoundary::Float64 = μ0-0.5
		UpperBoundary::Float64 = μ0+0.5
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

	end

	if debug
		n = sum( GetKPopulation(Phase,Parameters,K,v,μ,β;OptimizeBZ) )/D
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
	RenormalizeBands::Bool=true,
	OptimizeBZ::Bool=true,
	μ0::Float64=0.0,
	debug::Bool=false
)::Tuple{Dict{String,Float64},Float64}

Returns: (HF dictionary, chemical potential).
"""
function PerformHFStep(
	Phase::String,						# Mean field phase
	Parameters::Dict{String,Float64},	# Model parameters t,U,V
	K::Matrix{Vector{Float64}},			# BZ grid
	v0::Dict{String,Float64},			# HF initializers
	n::Float64,							# Density
	β::Float64;							# Inverse temperature
	Syms::Vector{String}=["s"],			# Gap function symmetries
	RenormalizeBands::Bool=true,			# Conditional renormalization of t
	OptimizeBZ::Bool=true,				# Conditional optimization of BZ
	Δn::Float64=1e-3,					# Density tolerance
	μ0::Float64=0.0,						# Initial guess for μ calculation
	debug::Bool=false,					# Debug mode
)::Tuple{Dict{String,Float64},Float64}

	v = copy(v0)	
	LxLy = prod(size(K))	
	μ = FindRootμ(Phase,Parameters,K,v0,n,β;Δn,μ0,RenormalizeBands,OptimizeBZ,debug)
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
			wk = GetWeight(q; Sym="S-MBZ",OptimizeBZ) # Avoid computational redundance
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
			wk = GetWeight(q;OptimizeBZ) # Avoid computational redundance
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
					sk = Δk/Ek * tanh(β*Ek/2)
					ck = ξk/Ek * tanh(β*Ek/2)
				end
				Δs -= wk * sk
				ΔS += wk * sk * StructureFactor("S",k)
				Δd += wk * sk * StructureFactor("d",k)
				gS += wk * (1-ck)/2 * StructureFactor("S",k)
				gd += wk * (1-ck)/2 * StructureFactor("d",k)
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
	RenormalizeBands::Bool=true,
	OptimizeBZ::Bool=true
)::Tuple{Dict{String,Any}, Dict{String,Any}}

Returns: (results dictionary, performance dictionary).
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
	RenormalizeBands::Bool=true,			# Conditional renormalization of t
	OptimizeBZ::Bool=true				# Conditional optimizatio of BZ
)::Tuple{Dict{String,Any}, Dict{String,Any}}

	if verbose
		@info "Running HF convergence algorithm" Phase Syms Parameters L n β
		@info "Convergence parameters" p Δv Δn g
	end

	# Reciprocal space discretization (normalized to 1)
	# Kx::Vector{Float64} = [kx for kx in -1:2/L[1]:1]
	# popfirst!(Kx)
	# Ky::Vector{Float64} = [ky for ky in -1:2/L[2]:1]
	# popfirst!(Ky)
	# K::Matrix{Vector{Float64}} = [ [kx,ky] for kx in Kx, ky in Ky ]
	K, Kx, Ky = GetKGrid(L)

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
				RenormalizeBands,
				OptimizeBZ,
				Δn,μ0=μ,
				debug,
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
				for key in keys(v0)
					current = v[key]
					previous = v0[key]
					if record
						Record[key] = vcat(Record[key],g*current + (1-g)*previous)
					end
					v[key] = g*current + (1-g)*previous
				end

				if debug
					@info "Initializer and current step after mixing" v0 v
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
		fMFT = GetFreeEnergy(Phase,Parameters,K,v0,n,μ,β;RenormalizeBands,OptimizeBZ)

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
