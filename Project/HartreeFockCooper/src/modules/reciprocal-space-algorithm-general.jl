#!/usr/bin/julia
using LinearAlgebra
using Roots
using Random

@doc raw"""
function GetHoppingEnergy(
    k::Vector{Float64}
)::Float64

Returns: the hopping energy and wavevector (`k[1]`, `k[2]`).

`GetHoppingEnergy` takes as input `k` (coordinate in k-space). It computes the 
hopping energy for the 2D monatomic regular square lattice.
"""
function GetHoppingEnergy(
    k::Vector{Float64},					# [kx, ky]
)::Float64
    return -2 * sum( cos.(k) )
end

@doc raw"""
function GetHamiltonian(
    k::Vector{Float64},
    m::Dict{String,Float64}
)::Matrix{Float64}

Returns: the 2×2 Nambu-Bogoliubov hamiltonian at wavevector (`k[1]`, `k[2]`).

`GetHamiltonian` takes as input `k` (coordinate in k-space) and `m` (dictionary 
of real HF parameters). It computes the contribution at wavevector (`k[1]`, 
`k[2]`) to the many-body second-quantized hamiltonian. The hamiltonian has
dimension (2,2) because Nambu spinors are used.
"""
function GetHamiltonian(
    k::Vector{Float64},					# [kx, ky]
	m::Dict{String,Float64}				# HF parameters
)::Matrix{Complex{Float64}}

    hk = zeros(Complex{Float64},2,2)    # Empty hamiltonian
    
	# Unpack HF parameters vector
	isnan( m["s"] )  ? Us = 0.0  : Us = m["s"]
	isnan( m["s*"] ) ? Vs = 0.0  : Vs = m["s*"] * sum( cos.(k) )
	isnan( m["px"] ) ? Vpx = 0.0 : Vpx = m["px"] * 1im * sin(k[1])
	isnan( m["py"] ) ? Vpy = 0.0 : Vpy = m["py"] * 1im * sin(k[2])
	isnan( m["d"] )  ? Vd = 0.0  : Vd = m["d"] * sum( cos.(k) .* [1,-1] )
	Δk = (Vs + Vpx + Vpy + Vd) - Us

    # Diagonal elements
    εk = GetHoppingEnergy(k)
    hk[1,1] = εk
    hk[2,2] = -εk
    
    # Off-diagonal elements
    hk[2,1] = -Δk
    hk[1,2] = -conj( Δk )

    return hk
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
    ε::Float64,                 		# Single-particle energy
    μ::Float64,                 		# Chemical potential
    β::Float64,                 		# Inverse temperature
)::Float64

    if β==Inf                   # Zero temperature distribution
        if ε<=μ
            return 1
        elseif ε>μ
            return 0
        end
    elseif β<Inf                # Finite temperature distribution
        return 1 / ( exp(β * (ε-μ)) + 1 )
    end
end

@doc raw"""
function GetQuasiParticlesEnergies(
    K::Matrix{Vector{Float64}},
    m0::Dict{String,Float64}
)::Matrix{Float64}

Returns:

`GetQuasiParticlesEnergies` takes as input `K` (k-points in the BZ), `m0` 
(dictionary of real HF parameters), `μ0` (chemical potential). It computes a
matrix of energies for quasiparticles at all wavevectors in BZ.
"""
function GetQuasiParticlesEnergies(
    K::Matrix{Vector{Float64}},			# BZ grid
    m0::Dict{String,Float64},        	# HF parameters
)::Matrix{Float64}
    
    @error "Deprecated."
    
    Ek 	= zeros(size(K))
	Δk = zeros(Complex{Float64}, size(K))
    
    for (i,k) in enumerate(K)
    
		# Unpack HF parameters vector
		isnan( m0["s"] )  ? Us = 0.0  : Us = m0["s"]
		isnan( m0["s*"] ) ? Vs = 0.0  : Vs = m0["s*"] * sum( cos.(k) )
		isnan( m0["px"] ) ? Vpx = 0.0 : Vpx = m0["px"] * 1im * sin(k[1])
		isnan( m0["py"] ) ? Vpy = 0.0 : Vpy = m0["py"] * 1im * sin(k[2])
		isnan( m0["d"] )  ? Vd = 0.0  : Vd = m0["d"] * sum( cos.(k) .* [1,-1] )
		Δk = (Vs + Vpx + Vpy + Vd) - Us
    end

	# Quasiparticles energies function
    E(ε::Float64, Δ::Complex{Float64}) = sqrt( ε^2 + abs(Δ)^2 )
    Ek .= E.( GetHoppingEnergy.(K), Δ )

    return Ek
end

@doc raw"""
function GetSingleStatePopulation(
    K::Matrix{Vector{Float64}},
    m::Dict{String,Float64},
    μ0::Float64,
    β::Float64;
    debug::Bool=false
)::Matrix{Float64}

Returns: matrix of single-particle k-states populations.

`GetSingleStatePopulation` takes as input `K` (k-points in the BZ), `m` 
(dictionary of real HF parameters), `μ` (chemical potential) and `β` (inverse
temperature). It computes a matrix of occupation numbers per each couple of
momenta.
"""
function GetSingleStatePopulation(
    K::Matrix{Vector{Float64}},			# BZ grid
    m::Dict{String,Float64},         	# HF parameters
    μ::Float64,                 		# Chemical potential
    β::Float64;                 		# Inverse temperature
    debug::Bool=false
)::Matrix{Float64}

    Nk = zeros(size(K))

    for (i,k) in enumerate(K)
    	hk = GetHamiltonian(k,m)
    	hk[1,1] -= μ # Include chemical potential
        hk[2,2] += μ # Include chemical potential
        F = eigen(hk)
        W = F.vectors
        E = real.( F.values ) 
        Tmp = ( abs(W[1,2])^2 - abs(W[1,1])^2 ) * tanh( β*E[2]/2 )
        if E[2]==0.0 && β==Inf
        	Tmp = 0.0	# Escape
        end
        
        if debug && isnan( Tmp )
        	@info "NaN detected" E W FermiDirac.(E,μ,β)
        end
        
        Nk[i] = 1 - Tmp
    end

    return Nk

end

@doc raw"""
function FindRootμ(
    K::Matrix{Vector{Float64}}
    m::Dict{String,Float64},
    nt::Float64,
    β::Float64;
    Δn::Float64=0.0
)::Float64

Returns: a variational estimation for the chemical potential at density `nt`.

`FindRootμ` takes as input `K` (k-points in the BZ), `m` (dictionary of real HF
parameters), `nt` (target density) and `β` (inverse temperature). It computes 
the root value for the chemical potential assuming the matrix in reciprocal 
space to be the Nambu-Bogoliubov hamiltonian obtained via the HF parameters `m`.
The optimization is performed using the `Roots.jl` library.
"""
function FindRootμ(
    K::Matrix{Vector{Float64}},			# BZ grid
    m::Dict{String,Float64},         	# HF parameters
    nt::Float64,                		# Target density
    β::Float64;                 		# Inverse temperature
    Δn::Float64=0.0,            		# Density tolerance
    debug::Bool=false
)::Float64

	if nt < 0 || nt > 1
		@error "Invalid target density. Choose 0 ≤ nt ≤ 1." nt
		return
	end

    D = 2 * prod(size(K))
    # Define function to be minimized
    δn(μ::Float64) = sum( GetSingleStatePopulation(K,m,μ,β) )/D - nt

    LowerBoundary = 0.0
    UpperBoundary = LowerBoundary
    if δn(LowerBoundary) >= 0
        while δn(LowerBoundary) >= 0
            if debug
                @warn "Moving down lower boundary"
            end
            LowerBoundary -= 1.0
        end
        UpperBoundary = LowerBoundary + 1.0
    elseif δn(LowerBoundary) < 0
        while δn(UpperBoundary) < 0
            if debug
                @warn "Moving up upper boundary"
            end
            UpperBoundary += 1.0
        end
        LowerBoundary = UpperBoundary - 1.0        
    end
    μ = find_zero(δn, (LowerBoundary, UpperBoundary))

    if debug
        n = sum( GetSingleStatePopulation(K,m,μ,β) )/D
        @info "Optimal chemical potential and density:" μ n    
    end

    return μ
end

@doc raw"""
function PerformHFStep(
    K::Matrix{Vector{Float64}},
    m0::Dict{String,Float64},
    U::Float64,
    V::Float64,
    n::Float64,
    β::Float64;
    debug::Bool=false
)::Dict{String,Float64}

Returns: HF estimation for Cooper instability parameter based on input.

`PerformHFStep` takes as input `K` (k-points in the BZ) , `m0` (dictionary of 
real HF initializers), 'U` (local interaction), `V` (non-local interaction), `n`
(density) and `β` (inverse temperature). It performs a  single step of the
iterative Hartree-Fock (HF) analysis starting over a 2D square lattice whose
dimensions are given by `L`. The function estimates `m` using as input the
Nambu-Bogoliubov hamiltonian obtained via the function `GetHamiltonian` and the
optimal chemical potential for the density `n` for this hamiltonian.
"""
function PerformHFStep(
    K::Matrix{Vector{Float64}},			# BZ grid
    m0::Dict{String,Float64},			# HF initializers
    U::Float64,							# Local interaction
    V::Float64,							# Non-local interaction
    n::Float64,                 		# Density
    β::Float64;                 		# Inverse temperature
    debug::Bool=false
)::Dict{String,Float64}

	ϕ = zeros(Complex{Float64}, size(K))
	m = copy(m0)
	
	LxLy = prod(size(K))	
    μ = FindRootμ(K,m0,n,β;debug)

    for (i,k) in enumerate(K)
        hk = GetHamiltonian(k,m0)
        hk[1,1] -= μ # Include chemical potential
        hk[2,2] += μ # Include chemical potential
        F = eigen(hk)
        W = F.vectors
        Wd = adjoint(W)
        E = real.( F.values ) # Strangely explicit realization is needed (?)
        	ϕ[i] += sum( [W[l,1] * Wd[2,l] * FermiDirac(E[l],0.0,β) for l in 1:2] )
    end
    
    display( ϕ )
    
  	# Structure factors
    fs(k::Vector{Float64}) = sum( cos.(k) )				# s*-wave
    fpx(k::Vector{Float64}) = -1im * sin(k[1])			# py-wave
    fpy(k::Vector{Float64}) = -1im * sin(k[2])			# px-wave
    fd(k::Vector{Float64}) = sum( cos.(k) .* [1,-1] )	# d-wave
    
    # Initialize response
    Us, Vs, Vpx, Vpy, Vd = [0.0 for _ in 1:5]
    if !isnan(m0["s"])
	    Us = U * sum(ϕ) / ( 2*LxLy )
    end

	CheckV = copy(m0)
	pop!(CheckV, "s")
	if any( .!isnan.(values(CheckV)) )
		Vk = zeros(Complex{Float64}, size(K))
		for (i,k) in enumerate(K)
			Vk[i] = sum( fs.([k] .- K) .* ϕ )
		end
		Vk .*= 2 * V / prod(size(K))
		Vs = real( sum(Vk .* fs.(K)) / LxLy	)			# s*-wave
		Vpx = real( sum(Vk .* fpx.(K)) / LxLy )			# px-wave
		Vpy = real( sum(Vk .* fpy.(K)) / LxLy )			# py-wave
		Vd = real( sum(Vk .* fd.(K)) / LxLy )			# d-wave
	end	
	
	# Clear cache or save data
    isnan( m0["s"] )  ? Us = NaN  : m["s"] = real( Us )
	isnan( m0["s*"] ) ? Vs = NaN  : m["s*"] = real( Vs )
	isnan( m0["px"] ) ? Vpx = NaN : m["px"] = real( Vpx )
	isnan( m0["py"] ) ? Vpy = NaN : m["py"] = real( Vpy )
	isnan( m0["d"] )  ? Vd = NaN  : m["d"] = real( Vd )
    
    return m
end

@doc raw"""
function RunHFAlgorithm(
	Syms::Vector{String},
	U::Float64,
    V::Float64,
    L::Vector{Int64},
    n::Float64,
    β::Float64,
    p::Int64,
    Δm::Dict{String,Float64},
    Δn::Float64,
    g::Float64;
    verbose::Bool=false,
    debug::Bool=false
)::Tuple{Dict{String,Float64}, Vector{Float64}, Float64}

Returns: Hartree-Fock (HF) estimation for `m`, sequential relative convergence 
parameters `Q` and computational time `ΔT`.

`RunHFAlgorithm` takes as input `Syms` (vector of symmetries to be simulated,
formatted as Strings) `U` (local interaction), `V` (non-local interaction) where
ratio by `t` is intended, `L` (square lattice dimensions), `n` (density), `β`
(inverse temperature), `p` (maximum number of HF iterations), `Δm` (tolerance on
each order parameter) and `Δn` (tolerance on density in chemical  potential 
estimation), `g` (mixing parameter). It performs an iterative HF analysis over a
2D square lattice whose dimensions are given by `L`. All of the function's
positional arguments must be positive defined. The algorithm is iterative and 
runs at most `p` times using at each reiteration the result of the previous 
computation.
"""
function RunHFAlgorithm(
	Syms::Vector{String},		# Gap function symmetries
    U::Float64,					# Local interaction
    V::Float64,					# Non-local interaction
    L::Vector{Int64},           # [Lx, Ly]
    n::Float64,                 # Density
    β::Float64,                 # Inverse temperature
    p::Int64,                   # Number of iterations
    Δm::Dict{String,Float64},   # Tolerance on each order parameter
    Δn::Float64,                # Tolerance on density difference
    g::Float64;                 # Mixing parameter
    verbose::Bool=false,
    debug::Bool=false
)::Tuple{Dict{String,Float64}, Vector{Float64}, Float64}

	# Initialize order parameters
	m0::Dict{String,Float64} = Dict([])
	m0["s"] = NaN
	m0["s*"] = NaN
	m0["px"] = NaN
	m0["py"] = NaN
	m0["d"] = NaN

	AllSyms = ["s", "s*", "px", "py", "d"]
	if !all( Syms[i] in AllSyms for i in 1:length(Syms) )
	    @error "Invalid symmetries. Plase use s, s*, px, py, d."
	    return
	end
	
	for Sym in Syms
		m0[Sym] = rand()	# Random initializer
	end
	    
    if verbose
        @info "Running Hartree-Fock convergence algorithm" Syms U V L n β
        @info "Convergence parameters" p Δm Δn
    end

    Kx = [kx for kx in -pi:2*pi/L[1]:pi]
    pop!(Kx)
    Ky = [ky for ky in -pi:2*pi/L[2]:pi]
    pop!(Ky)
    
    K::Matrix{Vector{Float64}} = [ [kx,ky] for kx in Kx, ky in Ky]

    m = copy(m0) # Shallow copy of values
    Qs = [0.0 for _ in 1:5]
    i = 1

    ΔT = @elapsed begin
        while i<=p

            if debug
                printstyled("\n---Step $i---\n", color=:yellow)
            end

            m = copy( PerformHFStep(K,m0,U,V,n,β;debug) )

			pp = [v for v in values( m0[Sym] for Sym in Syms )]	# Previous
			cc = [v for v in values( m[Sym] for Sym in Syms )]	# Current
			tt = [v for v in values( Δm[Sym] for Sym in Syms )]	# Tolerances
            Qs = abs.( cc .- pp ) ./ tt

            if all(Qs .<= 1)
                
                if verbose
                    printstyled("\n---Converged at step $i---\n", color=:green)
                end
                i = p+1
                
            elseif any(Qs .> 1)
            	
            	# Evaluate: it could be better to run another iteration only for
            	# symmetries with low convergence quality factor.
            	for (s,Sym) in enumerate(Syms)
            		cSym = m[Sym]
            		pSym = m0[Sym]
            		m[Sym] = g*cSym + (1-g)*pSym
            	end
            	
            	if debug
		            @info "Initializer and current step" m0 m
		        end
                i += 1      
                
            end    
            m0 = copy(m)
        end
    end
    
    if verbose
        if all(Qs .<= 1)
            @info "Algorithm has converged" m Qs
        elseif any(Qs .> 1)
            @info "Algorithm has not converged" m Qs
        end
    end

    return m,Qs,ΔT
end
