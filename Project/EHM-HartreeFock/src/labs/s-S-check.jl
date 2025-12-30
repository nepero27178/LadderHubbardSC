using DataFrames
using DelimitedFiles
using Plots

# Includer
PROJECT_ROOT = @__DIR__
PROJECT_ROOT *= "/.."   # Up to the effective root
# include(PROJECT_ROOT * "/setup/graphic-setup.jl")

t = 1.0
U = 10.0
V = 2.0
β = 100.0
Lx = 128

# Reciprocal space discretization (normalized to 1)
Kx::Vector{Float64} = [kx for kx in -1:2/Lx:1]
popfirst!(Kx)
Ky::Vector{Float64} = [ky for ky in -1:2/Lx:1]
popfirst!(Ky)
K::Matrix{Vector{Float64}} = [ [kx,ky] for kx in Kx, ky in Ky ]

Δs0 = 1.0
ΔS0 = 1.0
μ = 1.058

xDomain = -5.0:0.5:5.0
yDomain = -5.0:0.5:5.0
zzs = zeros(length(yDomain),length(xDomain))
zzS = zeros(length(yDomain),length(xDomain))
Itns = length(xDomain) * length(yDomain)

FilePathOut = PROJECT_ROOT * "/labs/s-S-check/U=$(U)_V=$(V)_μ=$(μ).txt"
mkpath(dirname(FilePathOut))
Header = "" #"x;y;zs;zS\n"
write(FilePathOut, Header)

Itn = 0
for (n,Δs0) in enumerate(xDomain), (m,ΔS0) in enumerate(yDomain)
	global Itn += 1
	print("\e[2K\e[1GRun ($Itn/$Itns)")
	Δs = 0.0
	ΔS = 0.0
	for (j,k) in enumerate(K)
		cxcy = cos(k[1]) + cos(k[2])
		ξk = -2*t*cxcy - μ
		Δk = Δs0 + ΔS0*cxcy
		Ek = sqrt( ξk^2 + abs(Δk)^2 )

		# Gap factor
		sk::Float64 = 0.0
		if Ek!=0.0 # Otherwise add nothing
			sk = Δk/Ek * tanh(β*Ek/2)
		end
		Δs -= sk
		ΔS += sk * cxcy
	end
	Δs *= U / (2*Lx^2)
	ΔS *= V / (Lx^2)

	zzs[m,n] = Δs
	zzS[m,n] = ΔS

	open(FilePathOut, "a") do io
		writedlm(io, [Δs0 ΔS0 Δs ΔS], ';')
	end
end