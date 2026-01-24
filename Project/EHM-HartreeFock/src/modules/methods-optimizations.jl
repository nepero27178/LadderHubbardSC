@doc raw"""
function GetWeight(
	k::Vector{Float64};
	Sym::String="S",
	OptimizeBZ::Bool=true
)::Int64

Returns: symmetry-structured weight for vector k (expressed in pi units!) in BZ.
"""
function GetWeight(
	k::Vector{Float64};					# [kx, ky] in pi units
	Sym::String="S",						# Symmetry structure
	OptimizeBZ::Bool=true				# Option to disable optimization
)::Int64

	# Weights
	wk::Int64 = 0
	kx, ky = k

	if !OptimizeBZ
        wk = 1
	elseif Sym=="S-MBZ" # Half-sized Brillouin Zone

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

	elseif Sym=="S"

        if (kx==0 && ky==0) || (kx==1 && ky==1)
            wk = 1
        elseif kx >= 0 && ky > 0
			if kx < 1 && ky < 1
				# Bulk (four times)
				wk = 4
			elseif kx==1 || ky==1
				# Edge (two times due to nesting)
				wk = 2
            end
        end
	end

	return wk

end

@doc raw"""
function GetUc(
    t::Float64,
    L::Vector{Int64},
    δ::Float64,
    β::Float64;
    Method::String=\"DisSum\",
    RenormalizeBands::Bool=true,
	OptimizeBZ::Bool=true
)::Float64

[...]
"""
function GetUc(
    t::Float64,							# Hopping amplitude
    # V::Float64,							# Non-local attraction
    L::Vector{Int64},					# [Lx, Ly]
    δ::Float64,							# Doping
    β::Float64;							# Inverse temperature
    Method::String="DisSum",				# Computational method
    RenormalizeBands::Bool=true,			# Conditional renormalization of t
	OptimizeBZ::Bool=true				# Conditional BZ optimization
)::Float64

    # Reciprocal space discretization (normalized to 1)
    Kx::Vector{Float64} = [kx for kx in -1:2/L[1]:1]
    popfirst!(Kx)
    Ky::Vector{Float64} = [ky for ky in -1:2/L[2]:1]
    popfirst!(Ky)
    K::Matrix{Vector{Float64}} = [ [kx,ky] for kx in Kx, ky in Ky ]

	Parameters::Dict{String,Float64} = Dict("t" => t)
	Fakev::Dict{String,Float64} = Dict("F" => 4.)
	μ = FindRootμ("Free",Parameters,K,Fakev,0.5+δ,β;RenormalizeBands,OptimizeBZ)

    if !in(Method, ["NumInt", "DisSum"])
        @error "Wrong method. Acceptable: \"NumInt\", \"DisSum\"."

    # Numerical integration
    elseif Method=="NumInt"

        F(x,m) = Elliptic.K( sqrt(1-x^2) ) * tanh( β/2 * (4*t*x - m) ) / (x -m/(4*t))
        DomainDown = (-1,0)
        ProblemDown = IntegralProblem(F, DomainDown, μ)
        SolutionDown = solve(ProblemDown, HCubatureJL())

        DomainUp = (0,1)
        ProblemUp = IntegralProblem(F, DomainUp, μ)
        SolutionUp = solve(ProblemUp, HCubatureJL())
        Uc = (2*pi)^2 / (SolutionUp.u + SolutionDown.u)

    # Discrete sum
    elseif Method=="DisSum"

        u::Float64 = 0.0
        for k in K .* pi
            ξk = GetHoppingEnergy(t,k) - μ
            if ξk != 0.0
                u += tanh(β/2 * ξk) / ξk
            elseif ξk == 0.0 # Handle correctly the zero limit
                u += β/2
            end
        end
        Uc = 2*prod(size(K)) / u

    end

    return Uc

end

@doc raw"""
function GetOptimalg(
    U::Float64,
    Uc::Float64;
    Δ::Float64=0.1
)::Float64

[...]
"""
function GetOptimalg(
    U::Float64,                         # Local repulsion
    Uc::Float64;                        # Critical U
    ΔU::Float64=1.0 					# U tolerance
)::Float64

    u::Float64 = U/Uc
    gc::Float64 = 2/(u+1)
#     if gc > Δ
#         return gc-Δ
#     elseif gc <= Δ
#         @warn "Δ=$(Δ) too large, tolerance ignored. Consider reducing Δ or U."
#         return gc/2
#     end
    return gc * (1-gc/2 * ΔU/Uc)

end
