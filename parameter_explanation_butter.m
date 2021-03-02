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
result_butter = zeros(1,length(truncation_range));


% Store the parameter best performance for exp filter
L0_ep_butter = zeros(25,1);

% Special Fourier reconstruction 
s_Fourier_special = fftshift(ifft(fftshift(S)));

% This is the standard for each truncation
ep1 =  median(abs(abs(s)-abs(s_Fourier_special)));
result_butter_order = zeros(1,length(truncation_range));
result_butter_Wn = zeros(1,length(truncation_range));
order_save = zeros(1,50);
order_cr = zeros(1,50);
for threshold = truncation_range
    
    %   update the waitbar
    waitbar1 = waitbar(threshold/total_step);

    % Truncation
    S(1:threshold) = zeros(1,threshold);
    S(N-threshold+1:N) = zeros(1,threshold);
    
    
    
    
    % Grid search both frame length and order
%     for alpha = 1:1:151
    for Wn = 0.01:0.01:0.99
        order_save = zeros(1,50);
        for order = 1:1:20
            
            reconstructionButterW = MyButter(order,Wn,S);  
            error_ff_butter = abs(reconstructionButterW) - abs(s);
            points_butter = find(abs(error_ff_butter) < ep1);
            % Save the performance of all order
            order_save(order) = length(points_butter);
        end
        % Find best order
        order_max = find(order_save == max(order_save),1 );
% %         % Use best order to reconstruct
        
        reconstructionButterW =MyButter(order_max,Wn,S) ;  
        error_ff_butter = abs(reconstructionButterW) - abs(s);
        points_butter = find(abs(error_ff_butter) < ep1);
        % Save the performance of all framelen
        L0_ep_butter(round(100*Wn)) = length(points_butter);
        % Save the correspondent power
        order_cr(round(100*Wn)) = order_max;
    end
    re_butterW = find(L0_ep_butter == max(L0_ep_butter),1 );
    
    reconstructionButterW = MyButter(order_cr(re_butterW),re_butterW/100,S) ; 
    result_butter_order(1,total_step+1-threshold) = order_cr(re_butterW);
    result_butter_Wn(1,total_step+1-threshold) = re_butterW/100;
    
    error_ff_butter = abs(reconstructionButterW) - abs(s);
    points_butter = find(abs(error_ff_butter) < ep1);
    % Save the performance of all framelen
    result_butter(1,total_step+1-threshold) = length(points_butter);
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
    order = 2;
    Wn = 0.28;
    
    reconstructionButterW = MyButter(order,Wn,S) ; 
    error_ff_butter = abs(reconstructionButterW) - abs(s);
    points_butter = find(abs(error_ff_butter) < ep1);
    % Save the performance of all framelen
    result_butter(2,total_step+1-threshold) = length(points_butter);
end
delete(waitbar1)



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
  
    reconstructionButterW = MyButter(4,0.14,S); 
    error_ff_butter = abs(reconstructionButterW) - abs(s);
    points_butter = find(abs(error_ff_butter) < ep1);
    % Save the performance of all framelen
    result_butter(3,total_step+1-threshold) = length(points_butter);
end
delete(waitbar1)
%%
plot(truncation_range/20,result_butter,'Linewidth',linew_zi)
xlabel('The Cut Off Frequency','FontSize',FontS)
ylabel('Ep0 metric value','FontSize',FontS)
legend('Grid search order and cut off frequency',' order = 2, Wn = 0.28;','order = 4  Wn = 0.14','FontSize',FontS-10)
%%
figure(2)
plot(truncation_range/20, result_butter_order,truncation_range/20, result_butter_Wn)

filename = 'butter_f.mat';
save(filename)
