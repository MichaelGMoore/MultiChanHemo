function h5put(outfile,datasetName,variable)
varname = inputname(3);
h5create(outfile,[datasetName,'/',varname],size(variable),'DataType',class(variable));
h5write(outfile,[datasetName,'/',varname],variable)

end

