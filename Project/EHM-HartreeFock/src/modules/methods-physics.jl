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
	μ::Float64,
	β::Float64;
	debug::Bool=false,
	RenormalizeBands::Bool=true
)::Matrix{Float64}

Returns: matrix of single-particle k-states populations.

`GetKPopulation` takes as input `Phase` (string specifying the mean-field phase,
the allowed are \"AF\", \"FakeAF\", \"SU-Singlet\", \"FakeSU-Singlet\",
\"SU-Triplet\", \"FakeSU-Triplet\"), `Parameters`  (dictionary of model
parameters containing `t`, `U`, `V`), `K` (k-points in the BZ), `v` (dictionary
of real HF parameters), `μ` (chemical potential) and `β` (inverse temperature).
It computes a matrix of occupation numbers per each couple of momenta. The
boolean option `RenormalizeBands' allows for choosing to renormalize or not
the hopping parameter.
"""
function GetKPopulation(
	Phase::String,						# Mean field phase
	Parameters::Dict{String,Float64},	# Model parameters t,U,V
	K::Matrix{Vector{Float64}},			# BZ grid
	v::Dict{String,Float64},				# HF parameters
	μ::Float64,							# Chemical potential
	β::Float64;							# Inverse temperature
	debug::Bool=false,
	RenormalizeBands::Bool=true,			# Conditional renormalization of t
	OptimizeBZ::Bool=true				# Conditional BZ optimization
)::Matrix{Float64}

	Sym = "S"
	if in(Phase, ["AF", "FakeAF"])
		Sym *= "-MBZ"
	end

	Nk = zeros(size(K))
	wk::Int64 = 0
	Ek::Float64 = 0.0

	# Compute free energy hopping shift
	# LxLy::Int64 = prod(size(K))
	# cc::Matrix{Float64} = StructureFactor.("S",K.*pi)
	# εε::Matrix{Float64} = GetHoppingEnergy.(Parameters["t"],K.*pi)
	# ff::Matrix{Float64} = FermiDirac.(εε,μ,β)
	# w::Float64 = sum(cc.*ff)/LxLy # Bare bands correction

	for (i,q) in enumerate(K)

		wk = GetWeight(q; Sym, OptimizeBZ) # Avoid computational redundance
		k = q .* pi # Important: multiply k by pi

		if Phase=="Free"

			t = Parameters["t"]
			# if RenormalizeBands
				# Conditional renormalization of bands
			# 	t -= w * Parameters["V"]
			# end
			εk = GetHoppingEnergy(t,k)
			Nk[i] = 2*FermiDirac(εk,μ,β)

		elseif in(Phase, ["AF", "FakeAF"]) && wk >= 1

			t = Parameters["t"]
			if RenormalizeBands
				# Conditional renormalization of bands
				t -= v["w0"] * Parameters["V"]
			end

			# Renormalized bands
			εk = GetHoppingEnergy(t,k)

			# Renormalized gap
			reΔk::Float64 = v["m"] * (Parameters["U"] + 8*Parameters["V"])
			imΔk::Float64 = 2*v["wp"]*Parameters["V"] * StructureFactor("S",k)

			# Renormalized gapped bands
			Ek = sqrt( εk^2 + reΔk^2 + imΔk^2 )

			# Local population
			Nk[i] = wk * (FermiDirac(-Ek,μ,β) + FermiDirac(Ek,μ,β))

		elseif in(Phase, ["SU-Singlet", "FakeSU-Singlet"]) && in(wk,[1,2,4])

			t = Parameters["t"]
			if RenormalizeBands && "gS" in keys(v)
				# Conditional renormalization of bands
				t -= v["gS"]/2 * Parameters["V"]
			end

			AllSyms = ["Δs", "ΔS", "Δd"]
			Fakev = copy(v)
			delete!(Fakev, "gS")
			delete!(Fakev, "gd")
			if !issubset(collect(keys(Fakev)), AllSyms)
				@error "Invalid set of symmetries. Please choose from $(AllSyms)."
			end

			# Free bands
			ξk::Float64 = GetHoppingEnergy(t,k) - μ
			if RenormalizeBands && "gd" in keys(v)
				ξk += Parameters["V"] * v["gd"] * StructureFactor("d",k)
			end

			# Gap
			Δk::Complex{Float64} = 0.0 + 1im * 0.0
			for (key, value) in v
				if !in(key, ["gS", "gd"])
					key = String(key)
					key = String(chop(key, head=1, tail=0))
					Δk += value * StructureFactor(key,k)
				end
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
function GetFreeEnergy(
	Phase::String,
	Parameters::Dict{String,Float64},
	K::Matrix{Vector{Float64}},
	v::Dict{String,Float64},
	μ::Float64,
	β::Float64;
	debug::Bool=false,
	RenormalizeBands::Bool=true
)Float64

Returns: Free energy density for the given Phase.

`GetFreeEnergy` takes as input `Phase` (string specifying the mean-field phase,
the allowed are \"AF\", \"FakeAF\", \"SU-Singlet\", \"FakeSU-Singlet\",
\"SU-Triplet\", \"FakeSU-Triplet\"), `Parameters`  (dictionary of model
parameters containing `t`, `U`, `V`), `K` (k-points in the BZ), `v` (dictionary
of real HF parameters), `μ` (chemical potential) and `β` (inverse temperature).
It sums the free energy contribution per each couple of momenta. The boolean
option `RenormalizeBands' allows for choosing to renormalize or not the hopping
parameter.
"""
function GetFreeEnergy(
	Phase::String,						# Mean field phase
	Parameters::Dict{String,Float64},	# Model parameters t,U,V
	K::Matrix{Vector{Float64}},			# BZ grid
	v::Dict{String,Float64},				# HF parameters
	n::Float64,							# Density
	μ::Float64,							# Chemical potential
	β::Float64;							# Inverse temperature
	RenormalizeBands::Bool=true,			# Conditional renormalization of t
	OptimizeBZ::Bool=true				# Conditional optimization of BZ
)::Float64

	Sym = "S"
	if in(Phase, ["AF", "FakeAF"])
		Sym *= "-MBZ"
	end

	# f0MFT::Float64 = 0.0 # TODO Integrate free energy computation
	fMFT::Float64 = 0.0
	LxLy::Int64 = prod(size(K))

	# Compute free energy hopping shift
	# μ0::Float64 = FindRootμ("Free",Parameters,K,v,n,β;RenormalizeBands,OptimizeBZ)
	# cc::Matrix{Float64} = StructureFactor.("S",K.*pi)
	# εε::Matrix{Float64} = GetHoppingEnergy.(Parameters["t"],K.*pi)
	# ff::Matrix{Float64} = FermiDirac.(εε,μ0,β)
	# w::Float64 = sum(cc.*ff)/LxLy # Bare bands correction

	for (i,q) in enumerate(K)

		wk = GetWeight(q; Sym, OptimizeBZ) # Avoid computational redundance
		k = q .* pi # Important: multiply k by pi

		if Phase=="Free"

			@error "Under construction"
			# t = Parameters["t"]
			# if RenormalizeBands
				# Conditional renormalization of bands
			# 	t -= w * Parameters["V"]
			# end
			# εk = GetHoppingEnergy(t,k)
			# f0MFT += 2/β * log( 1-FermiDirac(εk,μ0,β) )

		elseif in(Phase, ["AF", "FakeAF"]) && wk >= 1
			@error "Under construction"

		elseif in(Phase, ["SU-Singlet", "FakeSU-Singlet"]) && in(wk,[1,2,4])

			t = Parameters["t"]
			if RenormalizeBands && "gS" in keys(v)
				# Conditional renormalization of bands
				t -= v["gS"]/2 * Parameters["V"]
			end

			AllSyms = ["Δs", "ΔS", "Δd"]
			Fakev = copy(v)
			delete!(Fakev, "gS")
			delete!(Fakev, "gd")
			if !issubset(collect(keys(Fakev)), AllSyms)
				@error "Invalid set of symmetries. Please choose from $(AllSyms)."
			end

			# Free bands
			ξk::Float64 = GetHoppingEnergy(t,k) - μ
			if RenormalizeBands && "gd" in keys(v)
				ξk += Parameters["V"] * v["gd"] * StructureFactor("d",k)
			end

			# Gap
			Δk::Complex{Float64} = 0.0 + 1im * 0.0
			for (key, value) in v
				if !in(key, ["gS", "gd"])
					key = String(key)
					key = String(chop(key, head=1, tail=0))
					Δk += value * StructureFactor(key,k)
				end
			end

			# Renormalized gapped bands
			Ek = sqrt( ξk^2 + abs(Δk)^2 )

			# Local population
			tk::Float64 = 0.0
			if Ek!=0.0
				tk = abs(Δk)^2/Ek * tanh(β*Ek/2)
			end
			fMFT += wk * (ξk - Ek + tk + 2/β * log( 1-FermiDirac(Ek,0.0,β) ))

		elseif Phase=="Su/Triplet"
			@error "Under construction"
		end
	end

	fMFT = fMFT/LxLy + 2*n*μ
	if Phase=="Free"
		fMFT = f0MFT

	elseif Phase=="SU-Singlet"
		Corr::Float64=0.0
		for gSym in ["gS", "gd"]
			try
				Corr += abs(v[gSym])^2
			catch
			end
		end

		fMFT -= 2 * Parameters["V"] * Corr # Lagrange correction
	end

	return fMFT
end
