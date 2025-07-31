using LinearAlgebra

function PerformHFStep(
    t::Float64,
    k::Vector{Float64},    
    m0::Float64,
    Up::Bool
)
    e0 = -2 * t * ( cos(k[1]) + cos(k[2]) )
    σ = 2*Up-1
    h = [e0 m0; m0 -e0]
    F = eigen(h)
    u = F.vectors[:,1]
    v = F.vectors[:,2]
    return u*u + v*v
end
