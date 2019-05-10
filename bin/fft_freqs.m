function [ freqs ] = fft_freqs(N,fs)
%FFT_FREQS returns the frequencies associated with fft components for a
%vector of length 'nodes' at sampling frequency 'fs'
% 'freqs' will be a N-by-1 array
% 'N' is number of sample points
% 'fs' is the sample frequency


if ~mod(N,2)
    % start with even case
    freqs = cat(2,0:(N/2-1),(-N/2:-1))*fs/N;
else
    % odd case
    freqs = cat(2,0:(N-1)/2,-(N-1)/2:-1)*fs/N;
end
freqs = transpose(freqs);
