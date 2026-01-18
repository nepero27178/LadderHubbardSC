function ImportData(
    FilePathIn::String
)::DataFrame

	# Read data and create DataFrame
	DataIn = open(FilePathIn) do io
		readdlm(FilePathIn, ';', comments=false, '\n')
	end
	DF = DataFrame(DataIn[2:end,:], DataIn[1,:])
    DF = identity.(DF) # Format columns type

	return DF
end
