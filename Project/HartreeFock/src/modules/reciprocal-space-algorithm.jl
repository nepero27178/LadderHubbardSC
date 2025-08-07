#!/usr/bin/julia
using LinearAlgebra
# using Optim
using Roots

@doc raw"""
function GetHoppingEnergy(
    k::Vector{Float64},
    t::Float64
)::Float64

Returns: the hopping energy at hopping amplitude `t` and wavevector (`k[1]`, 
`k[2]`).

`GetHoppingEnergy` takes as input `k` (coordinate in k-space) and `t` (hopping 
amplitude). It computes the hopping energy for the 2D monatomic regular 
square lattice.
"""
function GetHoppingEnergy(
    k::Vector{Float64},         # [kx, ky]
    t::Float64,                 # t/U
)::Float64
    return -2 * t * ( cos(k[1]) + cos(k[2]) )
end

@doc raw"""
function GetHamiltonian(
    k::Vector{Float64},
    t::Float64,
    m::Float64
)::Matrix{Float64}

Returns: the 2×2 Nambu-Bogoliubov hamiltonian at wavevector (`k[1]`, `k[2]`).

`GetHamiltonian` takes as input `k` (coordinate in k-space), `t` (hopping
amplitude) and `m` (magnetization). It computes the contribution at
wavevector (`k[1]`, `k[2]`) to the many-body second-quantized hamiltonian. The
hamiltonian has dimension (2,2) because Nambu spinors are used.
"""
function GetHamiltonian(
    k::Vector{Float64},         # [kx, ky]
    t::Float64,                 # t/U
    m::Float64,                 # Mean-field magnetization
)::Matrix{Float64}
    hk = zeros(2,2)             # Empty hamiltonian

    # Diagonal elements
    εk = GetHoppingEnergy(k,t)
    hk[1,1] = εk
    hk[2,2] = -εk
    
    # Off-diagonal elements
    hk[1,2] = -m
    hk[2,1] = -m

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
    ε::Float64,                 # Single-particle energy
    μ::Float64,                 # Chemical potential
    β::Float64,                 # Inverse temperature
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
    Kx::Vector{Float64},
    Ky::Vector{Float64},
    t::Float64,
    m0::Float64,
    μ0::Float64
)::Matrix{Float64}

Returns:

`GetQuasiParticlesEnergies` takes as input
"""
function GetQuasiParticlesEnergies(
    Kx::Vector{Float64},        # BZ x grid
    Ky::Vector{Float64},        # BZ y grid
    t::Float64,                 # t/U
    m0::Float64,                # Magnetization
    μ0::Float64,                # Chemical potential
)::Matrix{Float64}
    
    Ek = zeros(length(Kx), length(Ky))
    #E(ε, μ, m) = sqrt( (ε-μ)^2 + m^2 )
    E(ε, m) = sqrt( ε^2 + m^2 )


    for (i,kx) in enumerate(Kx), (j,ky) in enumerate(Ky)
        k = [kx,ky]
        
        # Limit to MBZ, leave the rest to zero
        if abs(kx)+abs(ky) <= pi
            εk = GetHoppingEnergy(k,t)
            #Ek[i,j] = E(εk, μ0, m0)
            Ek[i,j] = E(εk, m0)
        end
    end

    return Ek
end

@doc raw"""
function GetSingleStatePopulation(
    Kx::Vector{Float64},
    Ky::Vector{Float64},
    t::Float64,
    m::Float64,
    μ0::Float64,
    β::Float64,
)::Matrix{Float64}

Returns: matrix of single-particle k-states populations.

`GetSingleStatePopulation` takes as input `Kx` (k-points in the x direction) ,
`Ky` (k-points in the y direction), `t` (hopping amplitude), `m` 
(magnetization), `μ` (chemical potential) and `β` (inverse temperature). 
It computes a matrix of occupation numbers per each couple of momenta, limitedly
to the Magnetic Brillouin Zone (defined by the condition |`kx`|+|`ky`|<π).
"""
function GetSingleStatePopulation(
    Kx::Vector{Float64},        # BZ x grid
    Ky::Vector{Float64},        # BZ y grid
    t::Float64,                 # t/U
    m::Float64,                 # Magnetization
    μ::Float64,                 # Chemical potential
    β::Float64,                 # Inverse temperature
)::Matrix{Float64}

    Nk = zeros(length(Kx), length(Ky))
    Ek = GetQuasiParticlesEnergies(Kx,Ky,t,m,μ)

    for (i,kx) in enumerate(Kx), (j,ky) in enumerate(Ky)
        # Limit to MBZ, don't touch the rest
        if abs(kx)+abs(ky) <= pi
            # Consider spin DoF
            Nk[i,j] = 2 * (FermiDirac(Ek[i,j],μ,β) + FermiDirac(-Ek[i,j],μ,β))
        end
    end

    return Nk

end

@doc raw"""
function FindRootμ(
    Kx::Vector{Float64},
    Ky::Vector{Float64},
    t::Float64,
    m::Float64,
    nt::Float64,
    β::Float64;
    Δn::Float64=0.0
)::Float64

Returns: a variational estimation for the chemical potential at density `nt`.

`FindRootμ` takes as input `Kx` (k-points in the x direction) , `Ky` (k-points 
in the y direction), `t` (hopping amplitude), `m` (magnetization), `nt` (target
density) and `β` (inverse temperature). It computes the root value for the 
chemical potential assuming the matrix in reciprocal space to be the Nambu-
Bogoliubov hamiltonian obtained via the parameters `t` and `m`. The 
optimization is performed using the `Roots.jl` library.
"""
function FindRootμ(
    Kx::Vector{Float64},        # BZ x grid
    Ky::Vector{Float64},        # BZ y grid
    t::Float64,                 # t/U
    m::Float64,                 # Magnetization
    nt::Float64,                # Target density
    β::Float64;                 # Inverse temperature
    Δn::Float64=0.0,            # Tolerance in density estimation
    debug::Bool=false
)::Float64

    D = 2 * length(Kx) * length(Ky)
    δn(μ) = sum( GetSingleStatePopulation(Kx,Ky,t,m,μ,β) )/D - nt

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
        println("Optimal chemical potential: $(μ0)")    
        n = sum( GetSingleStatePopulation(Kx,Ky,t,m,μ,β) )/D
        println("Obtained density: $(n)")
    end

    return μ
end

@doc raw"""
function PerformHFStep(
    Kx::Vector{Float64},
    Ky::Vector{Float64},
    t::Float64,
    m0::Float64,
    n::Float64,
    β::Float64;
    debug::Bool=false
)::Float64

Returns: HF estimation for AF instability parameter based on input.

`PerformHFStep` takes as input `Kx` (k-points in the x direction) , `Ky` 
(k-points in the y direction), `t` (hopping amplitude), `m0` (previous 
magnetization), `n` (density) and `β` (inverse temperature). It performs a 
single step of the iterative Hartree-Fock (HF) analysis starting over a 2D 
square lattice whose dimensions are given by `L`. The function estimates `m` 
using as input the Nambu-Bogoliubov hamiltonian obtained via the function 
`GetHamiltonian` and the optimal chemical potential for the density `n` for this
hamiltonian.
"""
function PerformHFStep(
    Kx::Vector{Float64},        # BZ x grid
    Ky::Vector{Float64},        # BZ y grid
    t::Float64,                 # t/U
    m0::Float64,                # Mean-field magnetization
    n::Float64,                 # Density
    β::Float64;                 # Inverse temperature
    debug::Bool=false
)::Float64

    m = 0
    μ = FindRootμ(Kx,Ky,t,m0,n,β)
    if debug
        println("μ=$μ") # Debug
    end

    for kx in Kx, ky in Ky
        k = [kx,ky]
        hk = GetHamiltonian(k,t,m0)# - μ * I(2) # Correct..?
        F = eigen(hk)
        W = F.vectors
        E = F.values
        # In FermiDirac set the chemical potential to zero (already in E)
        m += sum( [W[l,1] * W[2,l] * FermiDirac(E[l],μ,β) for l in 1:2] )
    end
    
    m /= length(Kx) * length(Ky) # Divide by 2?
    return m
end

@doc raw"""
function RunHFAlgorithm(
    t::Float64,
    L::Vector{Int64},
    n::Float64,
    β::Float64,
    p::Int64,
    Δm::Float64,
    Δn::Float64;
    verbose::Bool=false,
    debug::Bool=false
)::Tuple{Float64, Float64, Float64}

Returns: Hartree-Fock (HF) estimation for `m`, sequential relative convergence 
parameter `Q` and computational time `ΔT`.

`RunHFAlgorithm` takes as input `t` (hopping amplitude), `L` (square lattice
dimensions), `n` (density), `β` (inverse temperature), `p` (maximum number of
HF iterations), `Δm` (tolerance on magnetization) and `Δn` (tolerance on density
in chemical potential estimation). It performs an iterative HF  analysis over a 
2D square lattice whose dimensions are given by `L`. All the  function's
positional arguments must be positive defined. The algorithm is iterative and 
runs at most `p` times using at each reiteration the result of the previous 
computation.
"""
function RunHFAlgorithm(
    t::Float64,                 # t/U
    L::Vector{Int64},           # [Lx, Ly]
    n::Float64,                 # Density
    β::Float64,                 # Inverse temperature
    p::Int64,                   # Number of iterations
    Δm::Float64,                # Tolerance on magnetization
    Δn::Float64;                # Tolerance on density difference
    verbose::Bool=false,
    debug::Bool=false
)::Tuple{Float64, Float64, Float64}

    #TODO Assert that every variable used is positive-defined
    
    if verbose
        @info "Running Hartree-Fock convergence algorithm" t L n β
        @info "Convergence parameters" p Δm Δn
    end

    Kx = [kx for kx in -pi:2*pi/L[1]:pi]
    pop!(Kx)
    Ky = [ky for ky in -pi:2*pi/L[2]:pi]
    pop!(Ky)

    m = 1-n                    # Initializer
    m0 = m    
    Q = 0.0
    i=1

    ΔT = @elapsed begin
        while i<=p

            if debug
                println("Step $i")
            end

            m = PerformHFStep(Kx,Ky,t,m0,n,β;debug)

            if debug
                println("m0=$(m0), m=$(m)\n")
            end

            Q = abs(m-m0) / Δm

            if Q <= 1
                if verbose
                    printstyled("Converged at step $i\n", color=:green)
                end
                i = p+1
            elseif Q > 1
                m0 = m
                i += 1      
            end    
        end
    end
    
    if verbose
        if Q <= 1
            @info "Algorithm has converged" m Q
        elseif Q > 1
            @info "Algorithm has not converged" m Q
        end
    end

    return m,Q,ΔT
end
