function obj = applyLowPassFilter(obj,params)
% apply a Low-pass filter
% Use FFT method
% soft-edge filter
% use autoregression signal extension to remove edge effects

assert(isfield(params,'filter'),'filter parameters must be initialized before attempting to filter data')
assert(isfield(params.filter,'f0'),'must specify filter edge frequency')
f0 = params.filter.f0;
if isfield(params.filter,'width')
    width = params.filter.width;
else
    width = 0;
end

N = size(obj.y{1},1);
D = size(obj.y{1},2);

% estimate number of extra frames needed 
P = ceil(10*obj.fs/f0); % P stands for Pad
times = cat(1,(-P:-1)'/obj.fs,obj.times,obj.times(end)+(1:P)'/obj.fs);
freqs = fft_freqs(N+2*P,obj.fs); % frequencies after padding

% create fft filter
erfFilter = .5*(1-erf((sqrt(pi)/(2*width))*(abs(freqs)-f0)));
hardFilter = double(abs(freqs) <= f0);

myLPF = hardFilter;

% fill extra frames with autoregression model
for c = 1:obj.numCh
    xL = flipud(obj.y{c}(1:P,:));
    xR = obj.y{c}((end-P+1):end,:);
    aL = arburg(xL, P/2);
    aR = arburg(xR, P/2);
    yL = zeros(P,D);
    yR = zeros(P,D);
    for d = 1:D
    [~, zfL] = filter(-[0 aL(d,2:end)], 1, xL(1:P,d)); 
    [~, zfR] = filter(-[0 aR(d,2:end)], 1, xR(1:P,d)); 
    yL(1:P,d) = flipud(filter([0 0], -aL(d,:), zeros(1, P), zfL)); 
    yR(1:P,d) = filter([0 0], -aR(d,:), zeros(1, P), zfR); 
    end
    clear d
    temp = cat(1,yL,obj.y{c},yR); % this is the data with extrapolated frames
    temp = real(ifft(fft(temp).*myLPF)); % filtered
    obj.y{c} = temp((P+1):(end-P),:); % replace data with filtered data   
end
clear c

obj.filtered = [f0,width];

end

