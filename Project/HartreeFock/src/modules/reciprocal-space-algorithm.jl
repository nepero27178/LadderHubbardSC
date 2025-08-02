using LinearAlgebra
using Optim

@doc raw"""
function GetHoppingEnergy(
    k::Vector{Float64},
    t::Float64
)::Float64

Returns: the hopping energy at hopping amplitude `t` and wavevector (`k[1]`, 
`k[2]`).
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

Returns: the 2·2 Nambu-Bogoliubov hamiltonian for the square Hubbard lattice at
wavevector (`k[1]`, `k[2]`).

`GetHamiltonian` computes the contribution at wavevector (`k[1]`, `k[2]`) to the
many-body second-quantized hamiltonian. The hamiltonian has dimension (2,2)
because Nambu spinors are used.
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

Returns: the Fermi Dirac distribution at energy `ε`, chemical potential `μ` and
inverse temperature `β`.
"""
function FermiDirac(
    ε::Float64,                 # Single-particle energy
    μ::Float64,                 # Chemical potential
    β::Float64,                 # Inverse temperature
)::Float64
    return 1 / ( exp(β * (ε-μ)) + 1 )
end

@doc raw"""
function FindOptimalμ(
    Kx::Vector{Float64},
    Ky::Vector{Float64},
    t::Float64,
    m0::Float64,
    n0::Float64,
    β::Float64;
    Δn::Float64=0.0
)::Float64

Returns: a variational estimation for the chemical potential at density n0.

`FindOptimalμ` computes the optimal value for the chemical potential assuming
the matrix in reciprocal space to be the Nambu-Bogoliubov hamiltonian obtained 
via the parameters `t` and `m0`. The optimization is performed using the 
Optim.jl library.
"""
function FindOptimalμ(
    Kx::Vector{Float64},        # BZ x grid
    Ky::Vector{Float64},        # BZ y grid
    t::Float64,                 # t/U
    m0::Float64,                # Magnetization
    n0::Float64,                # Target density
    β::Float64;                 # Inverse temperature
    Δn::Float64=0.0             # Tolerance in density estimation
)::Float64

    Ek = zeros(length(Kx),length(Ky))
    for (i,kx) in enumerate(Kx), (j,ky) in enumerate(Ky)
        k = [kx,ky]
        Ek[i,j] = sqrt( GetHoppingEnergy(k,t)^2 + m0^2 )
    end
    δn(μ) = abs( sum(FermiDirac.(Ek,μ[1],β) + FermiDirac.(-Ek,μ[1],β)) - n0 )
    G = optimize(
        δn, [n0],
        method = Newton(),
        f_abstol = Δn
    )

    return G.minimum
end

@doc raw"""
function PerformHFStep(
    Kx::Vector{Float64},
    Ky::Vector{Float64},
    t::Float64,
    m0::Float64,
    n0::Float64,
    β::Float64
)::Float64

Returns: HF estimation for AF instability parameter based on input.

`PerformHFStep` performs a single step of the iterative Hartree-Fock (HF)
analysis starting over a 2D square lattice whose dimensions are given by `L`.
The function estimates `m` using as input the Nambu-Bogoliubov hamiltonian
obtained via the function `GetHamiltonian` and the optimal chemical potential
for the density `n0` for this hamiltonian.
"""
function PerformHFStep(
    Kx::Vector{Float64},        # BZ x grid
    Ky::Vector{Float64},        # BZ y grid
    t::Float64,                 # t/U
    m0::Float64,                # Mean-field magnetization
    n0::Float64,                # Density
    β::Float64,                 # Inverse temperature
)::Float64
    m = 0

    μ = FindOptimalμ(Kx,Ky,t,m0,n0,β)
    for kx in Kx, ky in Ky
        k = [kx,ky]
        hk = GetHamiltonian(k,t,m0) - μ * I(2)
        F = eigen(hk)
        W = F.vectors
        E = reverse(F.values)   # First +, then -
        m += sum( [W[l,1] * W[2,l] * FermiDirac(E[l],μ,β) for l in 1:2] )
    end
    
    m /= 2 * length(Kx) * length(Ky)
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
    verbose::Bool=false
)::Float64

Returns: HF converged estimation for AF instability parameter.

`RunHFAlgorithm` performs an iterative Hartree-Fock (HF) analysis over a 2D
square lattice whose dimensions are given by `L`. All the function's positional
arguments must be positive defined. The algorithm is iterative and runs at most
`p` times using at each reiteration the result of the previous computation.
"""
function RunHFAlgorithm(
    t::Float64,                 # t/U
    L::Vector{Int64},           # [Lx, Ly]
    n::Float64,                 # Density
    β::Float64,                 # Inverse temperature
    p::Int64,                   # Number of iterations
    Δm::Float64,                # Tolerance on magnetization
    Δn::Float64;                # Tolerance on density difference
    verbose::Bool=false
)::Float64

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
    i=1
    while i<=p
        m = PerformHFStep(Kx,Ky,t,m0,n,β)
        if abs(m-m0) <= Δm
            if verbose
                printstyled("Converged at step $i\n", color=:green)
            end
            i = p+1
        elseif abs(m-m0) > Δm
            m0 = m
            i += 1      
        end    
    end

    Q = abs(m-m0) / Δm
    
    if verbose
        if Q <= 1
            @info "Algorithm has converged" m Q
        elseif Q > 1
            @info "Algorithm has not converged" m Q
        end
    end

    return m
end
