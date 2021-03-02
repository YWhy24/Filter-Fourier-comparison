function beta = beta_calculation(signal,reconstruction)
N = length(signal);
% Initial
numerator = 0;
denominator = 0;
% Calculate Beta
signal_highpass = abs(laplace1d(signal));    
reconstruction_highpass = abs(laplace1d(reconstruction));
signal_highpass_mean = mean(signal_highpass);
reconstruction_highpass_mean = mean(reconstruction_highpass);
    
for i = 1:N 
    numerator = numerator + (signal_highpass(i)-signal_highpass_mean)*(reconstruction_highpass(i)-reconstruction_highpass_mean); 
end

for i = 1:N
    denominator =  denominator+ (signal_highpass(i)-signal_highpass_mean)^2 *(reconstruction_highpass(i)-reconstruction_highpass_mean)^2    ;   
end

denominator = sqrt(denominator);
beta = numerator/ denominator;
end

