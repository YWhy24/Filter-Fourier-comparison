% This is a code for filter comparisons using the ep0 metric.
% Author: Yue Wang
% Date: 19th Nov, 2020
% Exponential filter
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
result_exp = zeros(1,length(truncation_range));


% Store the parameter best performance for exp filter
L0_ep_exp = zeros(25,1);

% Special Fourier reconstruction
s_Fourier_special = fftshift(ifft(fftshift(S)));

% This is the standard for each truncation
ep1 =  median(abs(s-s_Fourier_special));

power_save = zeros(1,50);
power_cr = zeros(1,50);
for threshold = truncation_range
    
    %   update the waitbar
    waitbar1 = waitbar(threshold/total_step);
    
    % Truncation
    S(1:threshold) = zeros(1,threshold);
    S(N-threshold+1:N) = zeros(1,threshold);
    
    % Grid search both frame length and order
    for alpha = 1:1:151
        power_save = zeros(1,50);
        for power = 1:1:20
            exponentialfilter = zeros(1,N);
            index = 1;
            for i = -1:2/N:1-2/N      %Range from -1 to 1, step size match FFT matrix size
                exponentialfilter(index) = exp(-alpha*(sqrt(i^2))^power);
                index = index+1;
            end
            reconstruction_exp = fftshift(ifft(fftshift(S.*exponentialfilter)));
            error_ffre_exp = reconstruction_exp - s;
            points_exp = find(abs(error_ffre_exp) < ep1);
            % Save the performance of all order
            power_save(power) = length(points_exp);
        end
        % Find best order
        power_max = find(power_save == max(power_save),1 );
        % Use best order to reconstruct
        exponentialfilter = zeros(1,N);
        index = 1;
        for i = -1:2/N:1-2/N      %Range from -1 to 1, step size match FFT matrix size
            exponentialfilter(index) = exp(-alpha*(sqrt(i^2))^power_max);
            index = index+1;
        end
        reconstruction_exp = fftshift(ifft(fftshift(S.*exponentialfilter)));
        error_ffre_exp = reconstruction_exp - s;
        points_exp = find(abs(error_ffre_exp) < ep1);
        % Save the performance of all framelen
        L0_ep_exp(alpha) = length(points_exp);
        % Save the correspondent power
        power_cr(alpha) = power_max;
    end
    re_exp = find(L0_ep_exp == max(L0_ep_exp),1 );
    exponentialfilter = zeros(1,N);
    index = 1;
    for i = -1:2/N:1-2/N      %Range from -1 to 1, step size match FFT matrix size
        exponentialfilter(index) = exp(-re_exp*(sqrt(i^2))^power_cr(re_exp));
        index = index+1;
    end
    reconstruction_exp = fftshift(ifft(fftshift(S.*exponentialfilter)));
    error_ffre_exp = reconstruction_exp- s;
    points_exp = find(abs(error_ffre_exp) < ep1);
    % Save the performance of all framelen
    result_exp(1,total_step+1-threshold) = length(points_exp);
end

delete(waitbar1)
%%
% Define the Fourier Transform using math calculation
w = 10;
f = -w:2*w/N:w-2*w/N;
S = 80*sinc(4*f);

for threshold = truncation_range
    
    %   update the waitbar
    waitbar1 = waitbar(threshold/total_step);
    
    % Truncation
    S(1:threshold) = zeros(1,threshold);
    S(N-threshold+1:N) = zeros(1,threshold);

    exponentialfilter = zeros(1,N);
    index = 1;
    for i = -1:2/N:1-2/N      %Range from -1 to 1, step size match FFT matrix size
        exponentialfilter(index) = exp(-32*(sqrt(i^2))^8);
        index = index+1;
    end
    reconstruction_exp = fftshift(ifft(fftshift(S.*exponentialfilter)));
    error_ffre_exp = reconstruction_exp - s;
    points_exp = find(abs(error_ffre_exp) < ep1);
    % Save the performance of all framelen
    result_exp(2,total_step+1-threshold) = length(points_exp);
end


plot(truncation_range/20,result_exp,'Linewidth',linew_zi)
xlabel('The Cut Off Frequency','FontSize',FontS)
ylabel('Ep0 metric value','FontSize',FontS)



filename = 'exp_f.mat';
save(filename)
