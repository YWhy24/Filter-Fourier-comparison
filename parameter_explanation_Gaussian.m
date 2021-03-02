clear all

% Plotting related parameters
linew = 0.8;
linew_zi = 0.8;
FontS = 26;

% Signal Length Definition
width = 20;
fs =20;% sampling frequency
N = width*fs ;% 400
% Determining the x axis for the original signal
x = -width/2:width/N:width/2-width/N; %[-10:10]
% Define the signal
s = zeros(1, N);
t2 = (2*N/5)+1:(3*N/5);%[161:240]
s(t2) = 1;
% Define the Fourier Transform using math calculation
w = fs/2;
f = -w:2*w/N:w-2*w/N;
S = 4*fs*sinc(4*f);

% The total number of truncations tested
total_step = 199;

% the truncation range. 
truncation_range = 1:1:total_step;

% Store the result
result_Gau = zeros(3,length(truncation_range));

% Define sigma range for the Gaussian filter
number_of_sigma_tested = 150;
sigma_start = 1;
sigma_range = sigma_start:1:sigma_start+number_of_sigma_tested-1;
length_Gau = zeros(length(sigma_range),1);


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
    
    % Fourier reconstruction
    s_Fourier_re = fftshift(ifft(fftshift(S)));
  
       
    % Gaussian reconstruction    
    for sigma = sigma_range
        reconstruction_gaussian = G_filter(S,sigma);
        error_gau = abs(reconstruction_gaussian) - abs(s);
        points_Gau = find(abs(error_gau) < ep1);
        length_Gau(sigma,1) = length(points_Gau);
    end
    % Best sigma for Gaussian   
    Gaussian_best_sigma = median(find(length_Gau == max(length_Gau)));  
    reconstruction_gaussian_best = G_filter(S,Gaussian_best_sigma);
    error_gau = abs(reconstruction_gaussian_best) - abs(s);
    points_Gau_best = find(abs(error_gau) < ep1);
    result_Gau(1,N/2-threshold) = length(points_Gau_best); 
    
end
delete(waitbar1)

% Define the Fourier Transform using math calculation
w = 10;
f = -w:2*w/N:w-2*w/N;
S = 80*sinc(4*f);

% Start the loop to change threshold
for threshold = truncation_range
    %   update the waitbar
    waitbar1 = waitbar(threshold/total_step);

    % Truncation
    S(1:threshold) = zeros(1,threshold);
    S(N-threshold+1:N) = zeros(1,threshold);
    
    % Fourier reconstruction
    s_Fourier_re = fftshift(ifft(fftshift(S)));
  
   
    reconstruction_gaussian_best = G_filter(S,60);
    error_gau = abs(reconstruction_gaussian_best) - abs(s);
    points_Gau_best = find(abs(error_gau) < ep1);
    result_Gau(2,N/2-threshold) = length(points_Gau_best);
  
    
end
delete(waitbar1)

% Define the Fourier Transform using math calculation
w = 10;
f = -w:2*w/N:w-2*w/N;
S = 80*sinc(4*f);
% Start the loop to change threshold
for threshold = truncation_range
    %   update the waitbar
    waitbar1 = waitbar(threshold/total_step);

    % Truncation
    S(1:threshold) = zeros(1,threshold);
    S(N-threshold+1:N) = zeros(1,threshold);
    
    % Fourier reconstruction
    s_Fourier_re = fftshift(ifft(fftshift(S)));


    % Filters
   
    reconstruction_gaussian_best = G_filter(S,30);
    error_gau = abs(reconstruction_gaussian_best) - abs(s);
    points_Gau_best = find(abs(error_gau) < ep1);
    result_Gau(3,N/2-threshold) = length(points_Gau_best);
  
    
end
delete(waitbar1)

% Grid search
plot(truncation_range/20,result_Gau(1,:),'Linewidth',linew_zi)
hold on
% sigma = 60
plot(truncation_range/20,result_Gau(2,:),'Linewidth',linew_zi)
hold on
xlabel('The Cut Off Frequency','FontSize',FontS)
ylabel('Ep0 metric value','FontSize',FontS)
% sigma = 30
plot(truncation_range/20,result_Gau(3,:),'Linewidth',linew_zi)
hold off
title('The Gaussian Filter Parameter Comparison','FontSize',FontS )
legend('sigma Grid search','sigma = 60', 'sigma = 30', 'Location', 'best','FontSize',FontS)

filename = 'Gaussian_final';
save(filename)