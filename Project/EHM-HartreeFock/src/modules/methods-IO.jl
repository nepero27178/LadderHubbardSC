using DataFrames
using DelimitedFiles

function ImportData(
    FilePathIn::String
)::DataFrame

	# Read data and create DataFrame
	DataIn = open(FilePathIn) do io
		readdlm(FilePathIn, ';', comments=false, '\n')
	end
	DF = DataFrame(DataIn[2:end,:], DataIn[1,:])
    DF = identity.(DF) # Format columns type
    DF.Lx = Int64.(DF.Lx)

	return DF
end

	
function UnpackFilePath(
	FilePathIn::String
)::Tuple{String,String,Vector{String}}

	SubStrs = split(FilePathIn,'/')
	SetupStr = SubStrs[end-2]
	PhaseStr = SubStrs[end-1]
	SymsStr = SubStrs[end]
	
	Setup = split(SetupStr,'=')[2]
	Phase = split(PhaseStr,'=')[2]
	Syms = split(split(SymsStr,'=')[2],'.')[1]

	return String(Setup), String(Phase), [string(s) for s in Syms]
end
	
function ReshapeData(
	DF::DataFrame;
	xVar::String="U",
	yVar::String="V",
	zVar::String="v"
)::Tuple{Any, Any, Any}
	
	xx = unique(DF[!,xVar])
	yy = unique(DF[!,yVar])
	zz = reshape(DF[!,zVar],length(yy),length(xx))
	
	return xx, yy, zz
end

function EnlargeDF!(
	DF::DataFrame
)
	
	# List HF parameters
	ListHFPs::Vector{String} = [key for key in keys(
		eval(Meta.parse(
			DF.v[1]
		))
	)]
	vv = eval.(Meta.parse.(DF.v))
	QQ = eval.(Meta.parse.(DF.Q))
	for HF in ListHFPs
		DF[!,HF] = get.(vv,HF,"N/A")
		DF[!,"Q" * HF] = get.(QQ,HF,"N/A")
	end
	
	return select!(DF, Not(:v, :Q))
	
end