function Snew = getSessionAttributes(S)

if isfield(S,'session')
    S = rmfield(S,'session');
end
for s = 1:S.numSessions
    info = h5info(S.fileNames{s});
    S.session(s,1).animalID = h5readatt(S.fileNames{s},'/','animalID');
    S.session(s,1).sessionID = h5readatt(S.fileNames{s},'/','sessionID');
end

Snew = S;
end

