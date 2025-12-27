#!/usr/bin/julia
using LinearAlgebra
using Roots
using Random
LinearAlgebra.BLAS.set_num_threads(Threads.nthreads()) # Parallel optimization

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
	elseif Phase=="SU-Singlet"
		issubset(Syms, ["s", "S", "d"]) ? Singlet = true : throw(SymErr)
	elseif Phase=="SU-Triplet"
		issubset(Syms, ["px", "py", "p+", "p-"]) ? Triplet = true : throw(SymErr)
	end

	KeysList::Vector{String} = ["Empty"]
	AF ? KeysList = ["m", "w0", "wp"] : 0
	Singlet ? KeysList = ["Δ$(Sym)" for Sym in Syms] : 0
	Triplet ? KeysList = ["Δ$(Sym)" for Sym in Syms] : 0
	return KeysList

end

@doc raw"""
function GetRMPs(
	Phase::String,
	Syms::Vector{String}=["s"]
)::Vector{String}

Returns: Renormalized Model Parameters labels for the given Phase.
"""
function GetRMPs(
	Phase::String;						# Mean field phase
	Syms::Vector{String}=["s"]			# Gap function symmetries
)::Vector{String}

	AF = false
	Singlet = false
	Triplet = false
	SymErr = "Invalid symmetries. $(Syms) is incoherent with $(Phase)."
	if in(Phase, ["AF", "FakeAF"])
		AF = true
	elseif Phase=="SU-Singlet"
		issubset(Syms, ["s", "S", "d"]) ? Singlet = true : throw(SymErr)
	elseif Phase=="SU-Triplet"
		issubset(Syms, ["px", "py", "p+", "p-"]) ? Triplet = true : throw(SymErr)
	end

	KeysList::Vector{String} = ["Empty"]
	AF ? KeysList = ["reΔ_tilde", "imΔ_tilde", "t_tilde"] : 0
	Singlet ? KeysList = ["Δ$(Sym)" for Sym in Syms] : 0
	Triplet ? KeysList = ["Δ$(Sym)" for Sym in Syms] : 0
	return KeysList

end

@doc raw"""
function GetWeight(
	k::Vector{Float64};
	Sym::String="S"
)::Int64

Returns: symmetry-structured weight for vector k (expressed in pi units!) in BZ.
"""
function GetWeight(
	k::Vector{Float64};					# [kx, ky] in pi units
	Sym::String="S"						# Symmetry structure
)::Int64

	# Weights
	wk::Int64 = 0
	kx, ky = k

	if Sym=="S"
		if abs(kx)+abs(ky) <= 1
			if kx>=0 && ky>0 && kx+ky<1
				# Bulk (four times)
				wk = 4
			elseif kx+ky==1 && (kx!=1 && ky!=1)
				# Edge (two times due to nesting)
				wk = 2
			elseif kx==0 && (ky==0 || ky==1)
				# Special points
				wk = 1
			end
		end
	elseif Sym=="d"
		@warn "d-wave symmetry structure to be implemented."
		wk = 1
	end
	return wk
end

@doc raw"""
function StructureFactor(
	Sym::String,
	k::Vector{Float64}
)::Complex{Float64}

Returns: structure factor for symmetry `Sym` at wavevector (`k[1]`, `k[2]`).

`StructureFactor ' takes as input `Sym` (string specifying symmetry), whose
acceptable values are s, S, px, py, d, and `k` (coordinate in k-space). It
computes the relative symmetry structure facotor at the specified symmetry and
point in reciprocal space.
"""
function StructureFactor(
	Sym::String,                        # Symmetry
	k::Vector{Float64}                  # [kx, ky]
)::Complex{Float64}

	AllSyms = ["s", "S", "px", "py", "d"]
	if !in(Sym, AllSyms)
		@error "Invalid symmetries. Plase use s, S, px, py, d."
		return
	end

	if Sym=="s"
		return 1
	elseif Sym=="S"
		return sum( cos.(k) )           # s*-wave
	elseif Sym=="px"
		return 1im * sqrt(2) * sin(k[1])# py-wave
	elseif Sym=="py"
		return 1im * sqrt(2) * sin(k[2])# px-wave
	elseif Sym=="d"
		return sum( cos.(k) .* [1,-1] )	# d-wave
	end
end

@doc raw"""
function GetHoppingEnergy(
	t::Float64,
	k::Vector{Float64}
)::Float64

Returns: the hopping energy and wavevector (`k[1]`, `k[2]`).

`GetHoppingEnergy` takes as input `t` (hopping amplitude),  `k` (coordinate in 
k-space). It computes the  hopping energy for the 2D monatomic regular square 
lattice.
"""
function GetHoppingEnergy(
	t::Float64,							# Hopping amplitude
	k::Vector{Float64},					# [kx, ky]
)::Float64
	return -2 * t * StructureFactor("S",k)
end

@doc raw"""
function GetHamiltonian(
	Phase::String,
	Parameters::Dict{String,Float64},
	k::Vector{Float64},
	v::Dict{String,Float64};
	RenormalizeHopping::Bool=true
)::Matrix{Complex{Float64}}

Returns: the Nambu-Bogoliubov hamiltonian at wavevector (`k[1]`, `k[2]`).

`GetHamiltonian` takes as input `Phase` (string specifying the mean-field phase,
the allowed are \"AF\", \"FakeAF\", \"SU-Singlet\", \"SU-Triplet\"), `Parameters`
(dictionary of model parameters containing `t`, `U`, `V`), `k` (coordinate in 
k-space) and `v` (dictionary of real HF parameters). It computes the 
contribution at wavevector (`k[1]`, `k[2]`) to the many-body second-quantized 
hamiltonian. The hamiltonian dimension depends on the phase.

--- FINISH COMMENT HERE: HOW DIMENSION DEPENDS ON THE PHASE? ---
"""
function GetHamiltonian(
	Phase::String,						# Mean field phase
	Parameters::Dict{String,Float64},	# Model parameters t,U,V
	k::Vector{Float64},					# [kx, ky]
	v::Dict{String,Float64},				# HF parameters
	μ::Float64;							# Chemical potential
	RenormalizeHopping::Bool=true,		# Conditional renormalization of t
	# Syms::Vector{String}=["s"]			# Symmetry sector
)::Matrix{Complex{Float64}}

	if in(Phase, ["AF", "FakeAF"]) # Symmetry sector is physically s+s*
		
		# Empty hamiltonian
		hk = zeros(Complex{Float64},2,2)

		# Renormalized bands
		t = Parameters["t"]
		if RenormalizeHopping
			# Conditional renormalization of bands
			t -= v["w0"] * Parameters["V"]
		end
		εk = GetHoppingEnergy(t,k)
	
		# Renormalized gap		
		Δk = v["m"] * (Parameters["U"] + 8*Parameters["V"])
			+ 2 * 1im * v["wp"] * Parameters["V"] * StructureFactor("S",k)

		# Return matrix
		hk = [(εk-μ) -conj(Δk); -Δk (-εk-μ)]
		return hk

	elseif Phase=="SU-Singlet"

		# Empty hamiltonian
		hk = zeros(Complex{Float64},2,2)

		AllSyms = ["s", "S", "d"]
		if !issubset(keys(v), AllSyms)
			@error "Invalid set of symmetries. Please choose from $(AllSyms)."
		end

		# Free bands
		ξk = GetHoppingEnergy(t,k) - μ

		# Gap
		Δk = 0.0 + 1im * 0.0
		for key in keys(v)
			Δk += v[key] * StructureFactor(key,k)
		end

		# Return matrix
		hk = [ξk -conj(Δk); -Δk -ξk]
		return hk

	elseif Phase=="Su/Triplet"
		@error "Under construction"
		return
	end
end

@doc raw"""
function FermiDirac(
	ε::Float64,
	μ::Float64,
	β::Float64
)::Float64

Returns: the Fermi Dirac distribution.

`FermiDirac` takes as input `ε` (single particle energy), `μ` (chemical 
potential) and `β` (inverse temperature). It computes, at the specified point,
the Fermi-Dirac distribution 1/(exp(β(ε-μ))+1). For zero-temperature
simulations, set `β`=Inf to convert the distribution to a localized step
function.
"""
function FermiDirac(
	ε::Float64,						# Single-particle energy
	μ::Float64,						# Chemical potential
	β::Float64,						# Inverse temperature
)::Float64

	if β==Inf # Zero temperature distribution
		if ε<=μ
			return 1
		elseif ε>μ
			return 0
		end
	elseif β<Inf	 # Finite temperature distribution
		return 1 / ( exp(β * (ε-μ)) + 1 )
	end
end

@doc raw"""
function GetKPopulation(
	Phase::String,
	Parameters::Dict{String,Float64},
	K::Matrix{Vector{Float64}},
	v::Dict{String,Float64},
	μ0::Float64,
	β::Float64;
	debug::Bool=false,
	RenormalizeHopping::Bool=true
)::Matrix{Float64}

Returns: matrix of single-particle k-states populations.

`GetKPopulation` takes as input `Phase` (string specifying the mean-field phase,
the allowed are \"AF\", \"FakeAF\", \"SU-Singlet\", \"SU-Triplet\"),
`Parameters`  (dictionary of model parameters containing `t`, `U`, `V`), `K`
(k-points in the BZ), `v` (dictionary of real HF parameters), `μ` (chemical 
potential) and `β` (inverse temperature). It computes a matrix of occupation 
numbers per each couple of momenta. The boolean option `RenormalizeHopping' 
allows for choosing to renormalize or not the hopping parameter.
"""
function GetKPopulation(
	Phase::String,						# Mean field phase
	Parameters::Dict{String,Float64},	# Model parameters t,U,V
	K::Matrix{Vector{Float64}},			# BZ grid
	v::Dict{String,Float64},				# HF parameters
	μ::Float64,							# Chemical potential
	β::Float64;							# Inverse temperature
	debug::Bool=false,
	RenormalizeHopping::Bool=true       # Conditional renormalization of t
)::Matrix{Float64}

	Nk = zeros(size(K))
	t = Parameters["t"]
	if Phase=="FakeAF" && RenormalizeHopping
		# Conditional renormalization of bands
		t -= v["w0"] * Parameters["V"]
	end

	wk::Int64 = 0
	Ek::Float64 = 0.0
	for (i,q) in enumerate(K)
		wk = 1 # GetWeight(q) # Avoid computational redundance
		k = q .* pi # Important: multiply k by pi
		if in(Phase, ["AF", "FakeAF"]) && in(wk,[1,2,4])

			# Renormalized bands
			εk::Float64 = GetHoppingEnergy(t,k)

			# Renormalized gap
			reΔk::Float64 = v["m"] * (Parameters["U"] + 8*Parameters["V"])
			imΔk::Float64 = 2*v["wp"]*Parameters["V"] * StructureFactor("S",k)
			
			# Renormalized gapped bands
			Ek = sqrt( εk^2 + reΔk^2 + imΔk^2 )
			
			# Local population
			Nk[i] = wk * (FermiDirac(-Ek,μ,β) + FermiDirac(Ek,μ,β))

		elseif Phase=="SU-Singlet"

			AllSyms = ["Δs", "ΔS", "Δd"]
			if !issubset(keys(v), AllSyms)
				@error "Invalid set of symmetries. Please choose from $(AllSyms)."
			end

			# Free bands
			ξk::Float64 = GetHoppingEnergy(t,k) - μ

			# Gap
			Δk::Complex{Float64} = 0.0 + 1im * 0.0
			for (key, value) in v
				key = String(key)
				key = String(chop(key, head=1, tail=0))
				Δk += value * StructureFactor(key,k)
			end

			# Renormalized gapped bands
			Ek = sqrt( ξk^2 + abs(Δk)^2 )

			# Local population
			ck::Float64 = 0.0
			if Ek!=0.0
				ck = ξk/Ek
			end
			Nk[i] = wk * (1 - ck * tanh(β*Ek/2))

		elseif Phase=="Su/Triplet"
			@error "Under construction"
			return
		end
	end

	return Nk

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
	RenormalizeHopping::Bool=true
)::Float64

Returns: a variational estimation for the chemical potential at density `nt`.

`FindRootμ` takes as input `Phase` (string specifying the mean-field phase, the
allowed are \"AF\", \"FakeAF\", \"SU-Singlet\", \"SU-Triplet\"), `Parameters`
(dictionary of model parameters containing `t`, `U`, `V`), `K` (k-points in the 
BZ), `v` (dictionary of real HF parameters), `nt` (target density) and `β` 
(inverse temperature). It computes the root value for the chemical potential 
assuming the matrix in reciprocal space to be the Nambu-Bogoliubov hamiltonian 
obtained via the HF parameters `v`. The optimization is performed using the 
`Roots.jl` library. The boolean option `RenormalizeHopping' allows for choosing
to renormalize or not the hopping parameter.
"""
function FindRootμ(
	Phase::String,						# Mean field phase
	Parameters::Dict{String,Float64},	# Model parameters t,U,V
	K::Matrix{Vector{Float64}},			# BZ grid
	v::Dict{String,Float64},				# HF parameters
	nt::Float64,							# Target density
	β::Float64;							# Inverse temperature
	Δn::Float64=0.0,						# Density tolerance
	debug::Bool=false,
	RenormalizeHopping::Bool=true		# Conditional renormalization of t
)::Float64

	if nt < 0 || nt > 1
		@error "Invalid target density. Choose 0 ≤ nt ≤ 1." nt
		return
	end

	D::Int64 = 2 * prod(size(K))
	# Define function to be minimized
	δn(μ::Float64) = sum(
			GetKPopulation(Phase,Parameters,K,v,μ,β;RenormalizeHopping)
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
	RenormalizeHopping::Bool=true
)::Tuple{Dict{String,Float64},Float64}

Returns: HF estimation for Cooper instability parameter based on input.

`PerformHFStep` takes as input `Phase` (string specifying the mean-field phase, 
the allowed are \"AF\", \"FakeAF\", \"SU-Singlet\", \"SU-Triplet\"),
`Parameters` (dictionary of model parameters containing `t`, `U`, `V`), `K`
(k-points in the BZ), `v0` (dictionary of real HF initializers), `n` (density)
and `β` (inverse temperature). It performs a  single step of the iterative
Hartree-Fock (HF) analysis starting over a 2D square lattice whose dimensions
are given by `L`. Among the optional parameters is `Syms`, a strings vector
specifying the superconducting gap function symmetries to be simulated (allowed
are: \"s\", \"S\", \"d\" for the \"SU-Singlet\" phase, and \"px\", \"py\",
\"p+\", \"p-\" for the \"SU-Triplet\") phase. The function estimates `v` using
as input the Nambu-Bogoliubov hamiltonian obtained via the function
`GetHamiltonian` and the optimal chemical potential for the density `n` for this
hamiltonian. The boolean option `RenormalizeHopping' allows for choosing to
renormalize or not the hopping parameter.
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
	RenormalizeHopping::Bool=true		# Conditional renormalization of t
)::Tuple{Dict{String,Float64},Float64}

	v = copy(v0)	
	LxLy = prod(size(K))	
	μ = FindRootμ(Phase,Parameters,K,v0,n,β;debug)
	wk::Int64=0

	# Antiferromagnet
	if in(Phase, ["AF", "FakeAF"])
		m::Float64 = 0.0
		w0::Float64 = 0.0
		wpi::Float64 = 0.0
		
		t = Parameters["t"]
		if RenormalizeHopping
			# Conditional renormalization of bands
			t -= v["w0"] * Parameters["V"]
		end
		
		for (i,q) in enumerate(K)
			wk = GetWeight(q) # Avoid computational redundance
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

	elseif Phase=="SU-Singlet"
		#TODO Check: probably s*-wave wks can be used.

		Δs::Float64 = 0.0 # s-wave
		ΔS::Float64 = 0.0 # s*-wave
		Δd::Float64 = 0.0 # d-wave

		t = Parameters["t"]

		for (i,q) in enumerate(K)
			wk = 1 #TODO GetWeight(q) # Avoid computational redundance
			k = q .* pi # Important: multiply k by pi
			if in(wk,[1,2,4]) # Allowed weights
				# Free bands
				ξk = GetHoppingEnergy(t,k) - μ

				# Gap
				Δk::Float64 = 0.0
				for (key, value) in v
					key = String(key)
					key = String(chop(key, head=1, tail=0))
					Δk += value * StructureFactor(key,k)
				end

				# Renormalized gapped bands
				Ek::Float64 = sqrt( ξk^2 + abs(Δk)^2 )

				# Gap factor
				sk::Float64 = 0.0
				if Ek!=0.0 # Otherwise add nothing
					sk = wk * Δk/Ek * tanh(β*Ek/2)
				end
				Δs -= sk
				ΔS += sk * StructureFactor("S",k)
				Δd += sk * StructureFactor("d",k)
			end
		end

		Δs *= Parameters["U"] / (2*LxLy)
		ΔS *= Parameters["V"] / LxLy
		Δd *= Parameters["V"] / LxLy

		# Correct for sign fluctuations
		# Δs < 0 ? Δs = 0.0 : 0
		# ΔS < 0 ? ΔS = 0.0 : 0
		# Δd < 0 ? Δd = 0.0 : 0

		"Δs" in keys(v0) ? v["Δs"] = Δs : 0
		"ΔS" in keys(v0) ? v["ΔS"] = ΔS : 0
		"Δd" in keys(v0) ? v["Δd"] = Δd : 0

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
	RenormalizeHopping::Bool=true
)::Tuple{Dict{String,Float64}, Dict{String,Float64}, Float64, 
Dict{String,Vector{Float64}}}

Returns: Hartree-Fock (HF) estimation for `v`, sequential relative convergence 
parameters `Q` and computational time `ΔT`. If `record` is set to true, the
last output contains the entire evolution of the parameters.

`RunHFAlgorithm` takes as input `Phase` (string specifying the mean-field phase, 
the allowed are \"AF\", \"FakeAF\", \"SU-Singlet\", \"SU-Triplet\"), `Parameters`
(dictionary of model parameters containing `t`, `U`, `V`), `L` (square lattice 
dimensions), `n` (density), `β` (inverse temperature), `p` (maximum number of HF
iterations), `Δm` (tolerance on each order parameter) and `Δn` (tolerance on
density in chemical  potential estimation), `g` (mixing parameter). It performs
an iterative HF analysis over a 2D square lattice whose dimensions are given by
`L`. All of the function's positional arguments must be positive defined. The
algorithm is iterative and runs at most `p` times using at each reiteration the
result of the previous computation. The positional argument `v0` allows for
custom initialization of HF parameters. When left unspecified, `v0` elements
are randomly initialized. The `record` boolean positional parameter, when
activated, registers the entire evolution in the Vector{Float64} values of the
last output object. The boolean option `RenormalizeHopping' allows for choosing
to renormalize or not the hopping parameter.
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
	RenormalizeHopping::Bool=true		# Conditional renormalization of t
)::Tuple{Dict{String,Float64}, Dict{String,Float64}, Float64, Float64, Dict{String,Vector{Float64}}}

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
	elseif issubset(keys(v0i), HFPs) #all([ in(key,HFPs) for key in keys(v0i) ])
		for HFP in HFPs
			v0[HFP] = copy(v0i[HFP])
		end
	end
	v = copy(v0)    # Shallow copy of values
	Qs = copy(v0)   # Copy NaN keys

	# Initialize record matrix
	Record::Dict{String,Vector{Float64}} = Dict([
		key => [ v0[key] ] for key in HFPs
	])

	# Recursive run
	i = 1
	ΔT = @elapsed begin
		while i<=p

			if debug
				printstyled("\n---Step $i---\n", color=:yellow)
			end

			Results = PerformHFStep(
				Phase,
				Parameters,
				K,v0,n,β;
				Syms,
				debug,
				RenormalizeHopping
			)
			v = copy(Results[1])
			μ = Results[2]
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
				i = p+1

			elseif any([Qs[key] for key in keys(v0)] .>= 1)
				# g = rand() # Random mixing parameter
				for key in keys(v0)
					current = v[key]
					previous = v0[key]
					if record
						Record[key] = vcat(Record[key],current)
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

	if all([Qs[key] for key in keys(v0)] .<= 1)
		if verbose
			@info "Algorithm has converged." v Qs
		end
	elseif any([Qs[key] for key in keys(v0)] .> 1)
		if verbose
			@info "Algorithm has not converged - v saved as NaN." v Qs Phase
		end
		for key in keys(v0)
			v[key] = NaN
		end
	end

	return v,Qs,ΔT,μ,Record
end
