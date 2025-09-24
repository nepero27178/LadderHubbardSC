#!/usr/bin/julia
using LinearAlgebra
using Roots
using Random

@doc raw"""
function StructureFactor(
    Sym::String,
    k::Vector{Float64}
)::Complex{Float64}

Returns: structure factor for symmetry `Sym` at wavevector (`k[1]`, `k[2]`).

`StructureFactor ' takes as input `Sym` (string specifying symmetry), whose
acceptable values are s, s*, px, py, d, and `k` (coordinate in k-space). It
computes the relative symmetry structure facotor at the specified symmetry and
point in reciprocal space.
"""
function StructureFactor(
    Sym::String,                        # Symmetry
    k::Vector{Float64}                  # [kx, ky]
)::Complex{Float64}

    AllSyms = ["s", "s*", "px", "py", "d"]
	if !( Sym in AllSyms )
	    @error "Invalid symmetries. Plase use s, s*, px, py, d."
	    return
	end
    
    if Sym=="s"
        return 1
    elseif Sym=="s*"
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
    k::Vector{Float64}
)::Float64

Returns: the hopping energy and wavevector (`k[1]`, `k[2]`).

`GetHoppingEnergy` takes as input `k` (coordinate in k-space). It computes the 
hopping energy for the 2D monatomic regular square lattice.
"""
function GetHoppingEnergy(
    k::Vector{Float64},					# [kx, ky]
)::Float64
    return -2 * StructureFactor("s*",k)
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
    AllSyms = ["s", "s*", "px", "py", "d"]    
    Δk = Dict([
        Sym => 0.0+1im*0.0 for Sym in AllSyms
    ])

    for Sym in AllSyms
        # Sum only !NaN entries
        isnan(m[Sym]) ? 0 : Δk[Sym] = m[Sym] * StructureFactor(Sym,k)
    end
    Gap = sum( [v for v in values(Δk)] )

    # Diagonal elements
    εk = GetHoppingEnergy(k)
    hk[1,1] = εk
    hk[2,2] = -εk
    
    # Off-diagonal elements
    hk[2,1] = -Gap
    hk[1,2] = -conj( Gap )

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
    Syms::Vector{String},
    K::Matrix{Vector{Float64}},
    m0::Dict{String,Float64},
    U::Float64,
    V::Float64,
    n::Float64,
    β::Float64;
    debug::Bool=false
)::Dict{String,Float64}

Returns: HF estimation for Cooper instability parameter based on input.

`PerformHFStep` takes as input `Syms` (vector of symmetries to be simulated,
formatted as Strings), `K` (k-points in the BZ) , `m0` (dictionary of  real HF
initializers), 'U` (local interaction), `V` (non-local interaction), `n` 
(density) and `β` (inverse temperature). It performs a  single step of the
iterative Hartree-Fock (HF) analysis starting over a 2D square lattice whose
dimensions are given by `L`. The function estimates `m` using as input the
Nambu-Bogoliubov hamiltonian obtained via the function `GetHamiltonian` and the
optimal chemical potential for the density `n` for this hamiltonian.
"""
function PerformHFStep(
    Syms::Vector{String},		        # Gap function symmetries
    K::Matrix{Vector{Float64}},			# BZ grid
    m0::Dict{String,Float64},			# HF initializers
    U::Float64,							# Local interaction
    V::Float64,							# Non-local interaction
    n::Float64,                 		# Density
    β::Float64;                 		# Inverse temperature
    debug::Bool=false
)::Dict{String,Float64}

    AllSyms = ["s", "s*", "px", "py", "d"]    
    Factors = Dict([
        Sym => 0.0+1im*0.0 for Sym in AllSyms
    ])

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

#       # "By hands"		
#		ξk = hk[1,1]
#		Δk = -hk[2,1]
#		Ek = sqrt(ξk^2 + abs(Δk)^2)		
#		sθk = Δk/Ek
#		if Δk!=0.0
#			cζk = real(Δk)/Δk
#			sζk = imag(Δk)/Δk
#		elseif Δk==0.0
#			cζk = 0.0
#			sζk = 0.0
#		end
#		ϕ[i] = 1/2 * sθk * ( cζk + 1im * sζk ) * tanh(β * Ek/2)
    end

    for Sym in Syms
        isnan(m0[Sym]) ? 0 : Factors[Sym] = sum( StructureFactor.(Sym,K) .* ϕ )
        c = 1.0        
        if Sym=="s"
            c = -U/2
        elseif Sym in ["px", "py"]
            c = -1im * sqrt(2) * V
        elseif Sym in ["s*", "d"]
            c = V
        end        
        m[Sym] = real( Factors[Sym] * c/LxLy )
    end

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
)::Tuple{Dict{String,Float64}, Dict{String,Float64}, Float64}

    AllSyms = ["s", "s*", "px", "py", "d"]    
	if !all( Syms[i] in AllSyms for i in 1:length(Syms) )
	    @error "Invalid symmetries. Plase use s, s*, px, py, d."
	    return
	end
    m0 = Dict([
        Sym => NaN for Sym in AllSyms
    ])
	
    # Random initialization of a subsection of symmetries
	for Sym in Syms
		Sym=="s" ? c = -1 : c = +1
		m0[Sym] = c * rand()
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

    m = copy(m0)    # Shallow copy of values
    Qs = copy(m0)   # Copy NaN keys
    i = 1

    ΔT = @elapsed begin
        while i<=p

            if debug
                printstyled("\n---Step $i---\n", color=:yellow)
            end

            m = copy( PerformHFStep(Syms,K,m0,U,V,n,β;debug) )

            for Sym in Syms
            	cSym = copy( m[Sym] )
            	pSym = copy( m0[Sym] )
            	tSym = copy( Δm[Sym] )
                Qs[Sym] = abs( cSym-pSym ) / tSym 
            end

            if all([Qs[Sym] for Sym in Syms] .< 1)
                
                if verbose
                    printstyled("\n---Converged at step $i---\n", color=:green)
                end
                i = p+1
                
            elseif any([Qs[Sym] for Sym in Syms] .>= 1)
            	for Sym in Syms
            		cSym = copy( m[Sym] )
            		pSym = copy( m0[Sym] )
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
        if all([Qs[Sym] for Sym in Syms] .<= 1)
            @info "Algorithm has converged." m Qs
        elseif any([Qs[Sym] for Sym in Syms] .> 1)
            @info "Algorithm has not converged - m saved as NaN." m Qs Syms
            for Sym in Syms          
                m[Sym] = NaN
            end
        end
    end

    return m,Qs,ΔT
end
