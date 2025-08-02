using LinearAlgebra

function GetRealSpaceHamiltonian(
    t::Float64,         # t/U
    L::Vector{Int64},   # [Lx, Ly]
    m::Float64,         # Mean-field magnetization
)
    D = 2*prod(L)       # Hilbert space dimension
    H = zeros(D,D)      # Empty hamiltonian
    
    # Mean-field contributions
    s = +1              # Initializing sign factor
    for y in 1:L[2]
        for x in 1:L[1]
            j = x + (y-1)*L[1]
            H[2*j-1,2*j-1] = s * m
            H[2*j,2*j] = -s * m
            s *= -1        
        end
        s *= -1
    end

    # Kinetic contributions
    for α in 1:D
        Indices = [
            mod1(α-2*L[1],D),
            mod1(α-2,D),
            mod1(α+2,D),
            mod1(α+2*L[1],D),
        ]
        H[Indices,α] .= -t
    end

    return H
end

function PerformHFStep(
    t::Float64,         # t/U
    L::Vector{Int64},   # [Lx, Ly]
    m0::Float64,        # Mean-field magnetization
)
    D = 2*prod(L)       # Hilbert space dimension

    H = GetRealSpaceHamiltonian(t,L,m0)
    F = eigen(H)
    v = F.vectors

    # Half-filling code
    TmpSum = 0
    for l in 1:Int64(D/2)
        for λ in 1:Int64(D/2)
            Sign = (-1)^( (λ+1) - floor(Int64,λ/L[1]) * (L[1]-1) )
            TmpAdd = abs( v[2*λ-1,l] )^2 - abs( v[2*λ,l] )^2
            TmpSum += Sign * TmpAdd
        end
    end
    m = TmpSum/D
    return m
end

function PerformHFAlgorithm(

)

end
