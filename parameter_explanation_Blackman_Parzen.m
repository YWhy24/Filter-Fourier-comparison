% This is a code for filter comparisons using the ep0 metric.
% Author: Yue Wang
% Date: 19th Nov, 2020

clear all

% Waitbar definition
f2 = waitbar(0,'Processing');

% Plotting related parameters
linew = 0.8;
linew_zi = 0.8;
FontS = 26;

% Signal Length Definition
width = 20;
N = width*20;

% Determining the x axis for the original signal
x = -width/2:width/N:width/2-width/N;

% Define the signal
s = zeros(1, N);
t2 = (2*N/5)+1:(3*N/5);
s(t2) = 1;

% Define the Fourier Transform using math calculation
w = 10;
f = -w:2*w/N:w-2*w/N;
S = 80*sinc(4*f);

% The total number of truncations tested
total_step = 199;

% the truncation range. 
truncation_range = 1:1:total_step;

% Store the result
result_BP = zeros(2,length(truncation_range));

% Special Fourier reconstruction 
s_Fourier_special = fftshift(ifft(fftshift(S)));

% This is the standard for each truncation
ep1 =  median(abs(abs(s)-abs(s_Fourier_special)));

% Start the loop to change threshold
for threshold = truncation_range
    
    %   update the waitbar
    waitbar1 = waitbar(threshold/total_step);

    % Truncation
    S(1:threshold) = zeros(1,threshold);
    S(N-threshold+1:N) = zeros(1,threshold);


    % Blackman
    W_blackman = blackman(N,'symmetric');
    reconstruction_Blackman = fftshift(ifft(fftshift(W_blackman.*S')));
    error_Blackman = abs(reconstruction_Blackman') - abs(s);
    p1_Blackman = find(abs(error_Blackman) < ep1);
    result_BP(1,N/2-threshold) = length(p1_Blackman);
    
    % Parzen
    W_Parzen = parzenwin(N);
    reconstruction_Parzen = fftshift(ifft(fftshift(S.*W_Parzen')));%for 1D
    error_parzen = abs(reconstruction_Parzen) - abs(s);
    p1_parzen = find(abs(error_parzen) < ep1);
    result_BP(2,N/2-threshold) =length(p1_parzen);
    
end
delete(waitbar1)

plot(truncation_range/20,result_BP,'Linewidth',linew_zi)
xlabel('The Cut Off Frequency','FontSize',FontS)
ylabel('Ep0 metric value','FontSize',FontS)
% title('The Tukey parameter comparison','FontSize',FontS )
legend('Blackman', 'Parzen', 'Location', 'best','FontSize',FontS)



filename = 'BP_final';
save(filename)