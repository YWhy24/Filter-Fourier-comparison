function output = FFTandTruncation(signal,threshold)

N = length(signal);
S = fftshift(fft(fftshift(signal)));
S(1:threshold) = zeros(1,threshold);
S(N-threshold+1:N) = zeros(1,threshold);
output = S;

end

