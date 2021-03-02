clear all



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
result_tukey = zeros(2,length(truncation_range));

% Define sigma range for the Gaussian filter
number_of_sigma_tested = 150;
sigma_start = 1;
sigma_range = sigma_start:1:sigma_start+number_of_sigma_tested-1;
length_Gau = zeros(length(sigma_range),1);

ep0_value = zeros(length(truncation_range),1);

s_Fourier_special = fftshift(ifft(fftshift(S)));

ep1 =  median(abs(abs(s)-abs(s_Fourier_special)));

% Waitbar definition
waitbar1  = waitbar(0,'Processing r = 1');

% Start the loop to change threshold
for threshold = truncation_range
    %   update the waitbar
    waitbar1 = waitbar(threshold/total_step);

    % Truncation
    S(1:threshold) = zeros(1,threshold);
    S(N-threshold+1:N) = zeros(1,threshold);
    
    % Fourier reconstruction
    s_Fourier_re = fftshift(ifft(fftshift(S)));
    
    W_tukey = tukeywin(N,1);
    reconstruction_tukey = fftshift(ifft(fftshift(S.*W_tukey')));
    error_ffre_tukey = abs(reconstruction_tukey) - abs(s);
    points_tukey = find(abs(error_ffre_tukey) < ep1);
    
    % The first row is tukey
    result_tukey(1,N/2-threshold) =  length(points_tukey);
  
    
end
delete(waitbar1)


% Waitbar definition
waitbar2  = waitbar(0,'Processing r = 0.5');

% Define the Fourier Transform using math calculation
w = 10;
f = -w:2*w/N:w-2*w/N;
S = 80*sinc(4*f);
% Start the loop to change threshold
for threshold = truncation_range
    %   update the waitbar
    waitbar2 = waitbar(threshold/total_step);

    % Truncation
    S(1:threshold) = zeros(1,threshold);
    S(N-threshold+1:N) = zeros(1,threshold);
    
    % Fourier reconstruction
    s_Fourier_re = fftshift(ifft(fftshift(S)));

    W_tukey = tukeywin(N,0.5);
    reconstruction_tukey = fftshift(ifft(fftshift(S.*W_tukey')));
    error_ffre_tukey = abs(reconstruction_tukey) - abs(s);
    points_tukey = find(abs(error_ffre_tukey) < ep1);
    
    % The first row is tukey
    result_tukey(2,N/2-threshold) =  length(points_tukey);
      
end
delete(waitbar2)


plot(truncation_range/20,result_tukey,'Linewidth',linew_zi)
xlabel('The Cut Off Frequency','FontSize',FontS)
ylabel('Ep0 metric value','FontSize',FontS)
% title('The Tukey parameter comparison','FontSize',FontS )
legend('r = 1', 'r = 0.5', 'Location', 'best','FontSize',FontS)


filename = 'Tukey_final';
save(filename)
