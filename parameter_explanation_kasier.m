% This is a code for filter comparisons using the ep0 metric.
% Author: Yue Wang
% Date: 19th Nov, 2020

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
result_Kasier = zeros(3,length(truncation_range));


% Kasier filter parameter definition
number_of_Kasier_p_tested = 650;
Kasier_p_start = 1;
Kasier_p_range = Kasier_p_start:1:Kasier_p_start+number_of_Kasier_p_tested-1;
L0_ep_Kasier = zeros(number_of_Kasier_p_tested,1);
result_p = zeros(1,length(truncation_range));

% Store the parameter best performance for Savitzky Golay filter
L0_ep_sgol = zeros(25,1);

% Store the parameter best performance for Cheb2 filter
points_Cheby2 = zeros(25,1);

% Store the ep1 value for every truncation(The very small number for error comparison)
ep0_value = zeros(length(truncation_range),1);

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

    % Kasier Filter  
    for i = Kasier_p_range
        H_Kasier =  kaiser(N,i);
        reconstruction_Kasier = fftshift(ifft(fftshift(S.*H_Kasier')));
        error_Kasier = abs(abs(reconstruction_Kasier) - abs(s));
        p1_Kasier = find(abs(error_Kasier) < ep1);
        L0_ep_Kasier(i) = length(p1_Kasier);
    end
    re_Kasier = find(L0_ep_Kasier == max(L0_ep_Kasier),1 );
    H_Kasier =  kaiser(N,re_Kasier);
    reconstruction_Kasier = fftshift(ifft(fftshift(S.*H_Kasier')));
    error_Kasier = abs(abs(reconstruction_Kasier) - abs(s));
    p1_Kasier = find(abs(error_Kasier) < ep1);
    result_Kasier(1,N/2-threshold) =length(p1_Kasier);
    result_p(1,N/2-threshold) = re_Kasier;
end
delete(waitbar1)


% Define the Fourier Transform using math calculation
w = 10;
f = -w:2*w/N:w-2*w/N;
S = 80*sinc(4*f);
% Start the loop to change threshold
for threshold = truncation_range
    waitbar1 = waitbar(threshold/total_step);
    S(1:threshold) = zeros(1,threshold);
    S(N-threshold+1:N) = zeros(1,threshold);
    
    H_Kasier =  kaiser(N,50);
    reconstruction_Kasier = fftshift(ifft(fftshift(S.*H_Kasier')));
    error_Kasier =abs( abs(reconstruction_Kasier) - abs(s));
    p1_Kasier = find(abs(error_Kasier) < ep1);
    result_Kasier(2,N/2-threshold) =length(p1_Kasier);    
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
    S(1:threshold) = zeros(1,threshold);
    S(N-threshold+1:N) = zeros(1,threshold);

    H_Kasier =  kaiser(N,150);
    reconstruction_Kasier = fftshift(ifft(fftshift(S.*H_Kasier')));
    error_Kasier = abs(reconstruction_Kasier) - abs(s);
    p1_Kasier = find(abs(error_Kasier) < ep1);
    result_Kasier(3,N/2-threshold) =length(p1_Kasier);    
end
delete(waitbar1)

plot(truncation_range/20,result_Kasier,'Linewidth',linew_zi)
hold off
legend('Grid search beta','beta = 50','beta = 150','Location','best','FontSize',FontS-10)
xlabel('The Cut Off Frequency','FontSize',FontS)
ylabel('Ep0 metric value','FontSize',FontS)


filename = 'Kasier_final';
save(filename)