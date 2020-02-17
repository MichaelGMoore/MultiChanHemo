function spectra = readSpectra(fileName,numF,numR)
% Returns a structure containting the spectral data in HDF5 file 
spectra = struct;

info = h5info(fileName);
% find the group containing spectra
names = {info.Groups(strcmp({info.Groups.Name},'/spectra')).Datasets.Name};
% parse the names and read/organize the spectral data
for m = 1:length(names)
    data(m).type = names{m}(1);
    if strcmp(data(m).type,'F')
        temp = names{m}(2:(end-2));
        if isempty(temp)
            data(m).num = 1;
        else
            data(m).num = str2num(temp);
        end
        data(m).type =strcat(data(m).type,names{m}((end-1):end));
    elseif strcmp(data(m).type,'R')
        data(m).num = str2num(names{m}(end));    
    end
    data(m).raw = h5read(fileName,['/spectra/',names{m}]);
end
% find the min and max frequencies
fmin = Inf;
fmax = -Inf;
for m = 1:length(data)
    fmin = min(fmin,min(data(m).raw(:,1)));
    fmax = max(fmax,max(data(m).raw(:,1)));
end
% interpolate onto 1Hz grid
spectra.type = {data.type};
spectra.num = [data.num];
spectra.lambda(:,1) = floor(fmin):ceil(fmax);
for m = 1:length(data)
    spectra.data(:,m) = interp1(data(m).raw(:,1),data(m).raw(:,2),spectra.lambda,'pchip','extrap');
    % compute mean wavelength
    spectra.meanLambda(m) = sum(spectra.lambda.*spectra.data(:,m))/sum(spectra.data(:,m));
end



end

