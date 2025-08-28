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
    Δ::Float64
)::Matrix{Float64}

Returns: the 2×2 Nambu-Bogoliubov hamiltonian at wavevector (`k[1]`, `k[2]`).

`GetHamiltonian` takes as input `k` (coordinate in k-space) and `m` (d-wave real
HF parameter). It computes the contribution at wavevector (`k[1]`, `k[2]`) to 
the many-body second-quantized hamiltonian. The hamiltonian has  dimension (2,2)
because Nambu spinors are used.
"""
function GetHamiltonian(
    k::Vector{Float64},					# [kx, ky]
	Δs::Float64,
	Δse::Float64
)::Matrix{Float64}

    hk = zeros(Float64,2,2)
	Δk 	= Δs * Δse * sum( cos.(k) )

    # Diagonal elements
    εk = GetHoppingEnergy(k)
    hk[1,1] = εk
    hk[2,2] = -εk
    
    # Off-diagonal elements
    hk[2,1] = -Δk
    hk[1,2] = -Δk

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
    m::Float64,
    μ0::Float64,
    β::Float64;
    debug::Bool=false
)::Matrix{Float64}

Returns: matrix of single-particle k-states populations.

`GetSingleStatePopulation` takes as input `K` (k-points in the BZ), `m` (d-wave
real HF parameter), `μ` (chemical potential) and `β` (inverse temperature). It 
computes a matrix of occupation numbers per each couple of momenta.
"""
function GetSingleStatePopulation(
    K::Matrix{Vector{Float64}},			# BZ grid
    ms::Float64,         				# HF parameter
	mse::Float64,			         		# HF parameters
    μ::Float64,                 		# Chemical potential
    β::Float64;                 		# Inverse temperature
    debug::Bool=false
)::Matrix{Float64}

    Nk = zeros(size(K))

    for (i,k) in enumerate(K)
    	hk = GetHamiltonian(k,ms,mse)
    	hk[1,1] -= μ # Include chemical potential
        hk[2,2] += μ # Include chemical potential
        F = eigen(hk)
        W = F.vectors
        E = F.values
        
        if debug
        	@info "FermiDirac" E FermiDirac.(E,μ,β)
        end
        
        Tmp = 0.0 # Escape
        V = ( abs(W[1,1])^2 - abs(W[1,2])^2 ) * tanh( β*E[2]/2 )
        if E[2]==0.0 && β==Inf
        	V = Tmp # Escape
        end
        Nk[i] = 1 + V
    end

    return Nk

end

@doc raw"""
function FindRootμ(
    K::Matrix{Vector{Float64}}
    m::Float64,
    nt::Float64,
    β::Float64;
    Δn::Float64=0.0
)::Float64

Returns: a variational estimation for the chemical potential at density `nt`.

`FindRootμ` takes as input `K` (k-points in the BZ), `m` (d-wave real HF 
parameter), `nt` (target density) and `β` (inverse temperature). It computes the
root value for the chemical potential assuming the matrix in reciprocal space to
be the Nambu-Bogoliubov hamiltonian obtained via the HF parameters `m`. The
optimization is performed using the `Roots.jl` library.
"""
function FindRootμ(
    K::Matrix{Vector{Float64}},			# BZ grid
    ms::Float64,			         		# HF parameters
	mse::Float64,			         		# HF parameters
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
    δn(μ::Float64) = sum( GetSingleStatePopulation(K,ms, mse,μ,β) )/D - nt

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
        println("Optimal chemical potential: $(μ)")    
        n = sum( GetSingleStatePopulation(K,ms,mse,μ,β) )/D
        println("Obtained density: $(n)")
    end

    return μ
end

@doc raw"""
function PerformHFStep(
    K::Matrix{Vector{Float64}},
    m0::Float64,
    U::Float64,
    V::Float64,
    n::Float64,
    β::Float64;
    debug::Bool=false
)::Float64

Returns: HF estimation for Cooper instability parameter based on input.

`PerformHFStep` takes as input `K` (k-points in the BZ) , `m0` (d-wave real HF
initializer), 'U` (local interaction), `V` (non-local interaction), `n` 
(density) and `β` (inverse temperature). It performs a  single step of the
iterative Hartree-Fock (HF) analysis starting over a 2D square lattice whose
dimensions are given by `L`. The function estimates `m` using as input the
Nambu-Bogoliubov hamiltonian obtained via the function `GetHamiltonian` and the
optimal chemical potential for the density `n` for this hamiltonian.
"""
function PerformHFStep(
    K::Matrix{Vector{Float64}},			# BZ grid
    ms0::Float64,						# HF initializers
	mse0::Float64,						# HF initializers
    U::Float64,							# Local interaction
    V::Float64,							# Non-local interaction
    n::Float64,                 		# Density
    β::Float64;                 		# Inverse temperature
    debug::Bool=false
)::Tuple{Float64, Float64}

	# ϕ = zeros(Float64, size(K))
	ms = 0.0
	mse = 0.0
	LxLy = prod(size(K))
    μ = FindRootμ(K,ms0,mse0,n,β;debug)

    for (i,k) in enumerate(K)
        ek = GetHoppingEnergy(k)
        xk = ek - μ
        dk = ms0 + mse0 * sum( cos.(k) )
        Ek = sqrt(xk^2 + dk^2)
        	ϕ = 0.0        
        if Ek!==0.0
		    s2 = dk/Ek
		    c2 = xk/Ek
			ϕ = s2/2 * tanh(β*Ek/2)
		end
        ms += ϕ
		mse += ϕ * sum( cos.(k) )
    end    

    ms *= -U/ (2*LxLy)
	mse *= 2*V / LxLy
    return ms, mse
end

@doc raw"""
function RunHFAlgorithm(
	U::Float64,
    V::Float64,
    L::Vector{Int64},
    n::Float64,
    β::Float64,
    p::Int64,
    Δm::Float64,
    Δn::Float64,
    g::Float64;
    verbose::Bool=false,
    debug::Bool=false
)::Tuple{Vector{Float64}, Vector{Float64}, Float64}

Returns: Hartree-Fock (HF) estimation for `m`, sequential relative convergence 
parameters `Q` and computational time `ΔT`.

`RunHFAlgorithm` takes as input `U` (local interaction), `V` (non-local
interaction) where ratio by `t` is intended, `L` (square lattice dimensions),
`n` (density), `β` (inverse temperature), `p` (maximum number of HF iterations),
`Δm` (tolerance on each order parameter) and `Δn` (tolerance on density in
chemical  potential estimation), `g` (mixing parameter). It performs an 
iterative HF analysis over a  2D square lattice whose dimensions are given by
`L`. All of the function's positional arguments must be positive defined. The 
algorithm is iterative and  runs at most `p` times using at each reiteration the
result of the previous computation.
"""
function RunHFAlgorithm(
    U::Float64,					# Local interaction
    V::Float64,					# Non-local interaction
    L::Vector{Int64},           # [Lx, Ly]
    n::Float64,                 # Density
    β::Float64,                 # Inverse temperature
    p::Int64,                   # Number of iterations
    Δm::Float64,		        # Tolerance on each order parameter
    Δn::Float64,                # Tolerance on density difference
    g::Float64;                 # Mixing parameter
    verbose::Bool=false,
    debug::Bool=false
)::Tuple{Float64, Float64, Float64}

    #TODO Assert that every variable used is positive-defined
    
    if verbose
        @info "Running Hartree-Fock convergence algorithm" U V L n β
        @info "Convergence parameters" p Δm Δn
    end

    Kx = [kx for kx in -pi:2*pi/L[1]:pi]
    pop!(Kx)
    Ky = [ky for ky in -pi:2*pi/L[2]:pi]
    pop!(Ky)
    
    K::Matrix{Vector{Float64}} = [ [kx,ky] for kx in Kx, ky in Ky]

    ms = rand()
	mse = rand()
	ms0 = ms
	mse0 = mse    
    Q = 0.0
    i = 1

    ΔT = @elapsed begin
        while i<=p

            if debug
                println("Step $i")
            end

            ms, mse = PerformHFStep(K,ms0,mse0,U,V,n,β;debug)

            if debug
                println("ms0=$(ms0), mse0=$(mse0)")
                println("ms=$(ms), mse=$(mse)\n")
            end

            Qs = abs( ms-ms0 ) / Δm
            Qse = abs( mse-mse0 ) / Δm

            if all([Qs, Qse] .<= 1)
                if verbose
                    printstyled("Converged at step $i\n", color=:green)
                end
                i = p+1
            elseif any([Qs, Qse] .> 1)
                ms0 = g*ms + (1-g)*ms0
                mse0 = g*mse + (1-g)*mse0
                i += 1      
            end    
        end
    end
    
    if verbose
        if all([Qs, Qse] .<= 1)
            @info "Algorithm has converged" m Q
        elseif any([Qs, Qse] .> 1)
            @info "Algorithm has not converged" m Q
        end
    end

    return ms,mse,Q,ΔT
end
