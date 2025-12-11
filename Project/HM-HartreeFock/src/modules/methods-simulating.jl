#!/usr/bin/julia
using LinearAlgebra
using Roots
using Random
LinearAlgebra.BLAS.set_num_threads(Threads.nthreads()) # Parallel optimization

function GetHFPs(
    Phase::String
)::Vector{String}

    KeysList::Dict{String,Vector{String}} = Dict([
        "AF" => ["m"]
    ])

    return KeysList[Phase]
end

function GetRMPs(
    Phase::String
)::Vector{String}

    KeysList::Dict{String,Vector{String}} = Dict([
        "AF" => ["Δ"]
    ])

    return KeysList[Phase]
end

function GetWeight(
    k::Vector{Float64}                  # [kx, ky] in pi units
)::Int64
    wk::Int64 = 0
    kx, ky = k
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
    return wk
end

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

function GetHoppingEnergy(
	t::Float64,							# Hopping amplitude
    k::Vector{Float64},					# [kx, ky]
)::Float64
    return -2 * t * StructureFactor("s*",k)
end

"""
function GetHamiltonian(
	Phase::String,						# Mean field phase
	Parameters::Dict{String,Float64},	# Model parameters t,U,V
    k::Vector{Float64},					# [kx, ky]
	v::Dict{String,Float64};			# HF parameters
)::Matrix{Complex{Float64}}

	if in(Phase, ["AF"])
		
		# Empty hamiltonian
		hk = zeros(Complex{Float64},2,2)

		t = Parameters["t"]
		εk = GetHoppingEnergy(t,k)
	
		# Renormalized gap		
		Δk = v["m"] * Parameters["U"]
		
        # Return matrix
		hk = [εk -Δk; -Δk -εk]
		return hk
        
    end

end
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

function GetKPopulation(
	Phase::String,						# Mean field phase
	Parameters::Dict{String,Float64},	# Model parameters t,U,V
    K::Matrix{Vector{Float64}},			# BZ grid
    v::Dict{String,Float64},         	# HF parameters
    μ::Float64,                 		# Chemical potential
    β::Float64;                 		# Inverse temperature
    debug::Bool=false,
)::Matrix{Float64}

    Nk = zeros(size(K))
    t = Parameters["t"]
    
    wk::Int64=0
    for (i,q) in enumerate(K)
        wk = GetWeight(q) # Avoid computational redundance
        k = q .* pi # Important: multiply k by pi
		if in(Phase, ["AF", "FakeAF"]) && in(wk,[1,2,4])
			# Bands
			εk::Float64 = GetHoppingEnergy(t,k)
		    
		    # Gap		
		    Δk::Float64 = v["m"] * Parameters["U"]
			
			# Gapped bands
			Ek::Float64 = sqrt( εk^2 + Δk^2 )
			
			# Fermi-Dirac factor
			Nk[i] = 2 * wk * (FermiDirac(-Ek,μ,β) + FermiDirac(Ek,μ,β))
		end
    end

    return Nk

end

function FindRootμ(
	Phase::String,						# Mean field phase
	Parameters::Dict{String,Float64},	# Model parameters t,U,V
    K::Matrix{Vector{Float64}},			# BZ grid
    v::Dict{String,Float64},         	# HF parameters
    nt::Float64,                		# Target density
    β::Float64;                 		# Inverse temperature
    Δn::Float64=1e-4,            		# Density tolerance
    debug::Bool=false,
)::Float64

	if nt < 0 || nt > 1
		@error "Invalid target density. Choose 0 ≤ nt ≤ 1." nt
		return
	end

    D::Int64 = 2 * prod(size(K))
    # Define function to be minimized
    δn(μ::Float64) = sum( 
            GetKPopulation(Phase,Parameters,K,v,μ,β)
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

function PerformHFStep(
	Phase::String,						# Mean field phase
	Parameters::Dict{String,Float64},	# Model parameters t,U,V
    K::Matrix{Vector{Float64}},			# BZ grid
    v0::Dict{String,Float64},			# HF initializers
    n::Float64,                 		# Density
    β::Float64;                 		# Inverse temperature
	debug::Bool=false,					# Debug mode
)::Tuple{Dict{String,Float64},Float64}

	v = copy(v0)	
	LxLy = prod(size(K))	
    μ = FindRootμ(Phase,Parameters,K,v0,n,β;debug)

    # Antiferromagnet
	if in(Phase, ["AF", "FakeAF"])
		m::Float64 = 0.0
		t = Parameters["t"]
		
        wk::Int64=0
		for (i,q) in enumerate(K)
            wk = GetWeight(q) # Avoid computational redundance
            k = q .* pi # Important: multiply k by pi
            if in(wk,[1,2,4]) # Allowed weights
                # Bands
    			εk::Float64 = GetHoppingEnergy(t,k)
                
                # Gap		
                Δk::Float64 = v["m"] * Parameters["U"]
                
                # Renormalized gapped bands
                Ek::Float64 = sqrt( εk^2 + Δk^2 )
                
                # Fermi-Dirac factor
                FDk = FermiDirac(-Ek,μ,β) - FermiDirac(Ek,μ,β)
                
                if Ek!=0.0 # Otherwise add nothing
                    m += wk * Δk/Ek * FDk
                end
            end
            wk = 0
		end
		
		v["m"] = m/(2*LxLy)

	end

    return v, μ
end

function RunHFAlgorithm(
    Phase::String,						# Mean field phase
	Parameters::Dict{String,Float64},	# Model parameters t,U,V
    L::Vector{Int64},                   # [Lx, Ly]
    n::Float64,                         # Density
    β::Float64,                         # Inverse temperature
    p::Int64,                           # Number of iterations
    Δv::Dict{String,Float64},           # Tolerance on each order parameter
    Δn::Float64,                        # Tolerance on density difference
    g::Float64;                         # Mixing parameter
    v0i::Dict=Dict([]),                 # Initializers
    verbose::Bool=false,
    debug::Bool=false,
    record::Bool=false,
)::Tuple{Dict{String,Float64}, Dict{String,Float64}, Float64, Float64, Dict{String,Vector{Float64}}}

    if verbose
        @info "Running HF convergence algorithm" Phase Parameters L n β
        @info "Convergence parameters" p Δv Δn g
    end

    # Reciprocal space discretization (normalized to 1)
    Kx::Vector{Float64} = [kx for kx in -1:2/L[1]:1]
    popfirst!(Kx)
    Ky::Vector{Float64} = [ky for ky in -1:2/L[2]:1]
    popfirst!(Ky)
    K::Matrix{Vector{Float64}} = [ [kx,ky] for kx in Kx, ky in Ky ]
    
    # Get Hartree Fock Parameters labels
    HFPs = GetHFPs(Phase)

    # Initialize HF dictionaries
    v0::Dict{String,Float64} = Dict([])
    μ::Float64 = 0.0

    if v0i==Dict([])
        for HFP in HFPs
            v0[HFP] = rand()
        end
    elseif all([ in(key,HFPs) for key in keys(v0i) ])
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
                debug,
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
    
    if verbose
        if all([Qs[key] for key in keys(v0)] .<= 1)
            @info "Algorithm has converged." v Qs
        elseif any([Qs[key] for key in keys(v0)] .> 1)
            @info "Algorithm has not converged - v saved as NaN." v Qs Phase
            for key in keys(v0)
                v[key] = NaN
            end
        end
    end

    return v,Qs,ΔT,μ,Record
end
