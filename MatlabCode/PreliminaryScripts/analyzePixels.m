S.fileNames = scanFolder(S.dataFolder);
S.numSessions = length(S.fileNames);
S = getSessionAttributes(S);
S = processPixels(S);
outName = ['S',S.dataSetName];
assignin('base',outName,S);
save(outName,outName)
clear outName