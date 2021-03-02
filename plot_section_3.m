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


% Store the parameter best performance for Savitzky Golay filter
L0_ep_sgol = zeros(25,1);

% Store the parameter best performance for Cheb2 filter
points_Cheby2 = zeros(25,1);

% Store the ep1 value for every truncation(The very small number for error comparison)
ep0_value = zeros(length(truncation_range),1);

% Special Fourier reconstruction
s_Fourier_special = fftshift(ifft(fftshift(S)));

% This is the standard for each truncation
ep1 =  median(abs(s-s_Fourier_special));
threshold = 120;


% Truncation
S(1:threshold) = zeros(1,threshold);
S(N-threshold+1:N) = zeros(1,threshold);


xf = -10:20/400:10-20/400;
N= 400;


W = hamming(400);
plot(xf,W,'r','LineWidth',1.3)
hold on

plot(xf,parzenwin(N),'-ro','LineWidth',1.3,'MarkerIndices',1:7:N)
hold on



for i = Kasier_p_range
    H_Kasier =  kaiser(N,i);
    reconstruction_Kasier = fftshift(ifft(fftshift(S.*H_Kasier')));
    error_Kasier = reconstruction_Kasier - s;
    p1_Kasier = find(abs(error_Kasier) < ep1);
    L0_ep_Kasier(i) = length(p1_Kasier);
end
re_Kasier = find(L0_ep_Kasier == max(L0_ep_Kasier),1 );
H_Kasier =  kaiser(N,re_Kasier);
reconstruction_Kasier = fftshift(ifft(fftshift(S.*H_Kasier')));
error_Kasier = abs(abs(reconstruction_Kasier) - abs(s));
p1_Kasier = find(abs(error_Kasier) < ep1);
length(p1_Kasier);

plot(xf,H_Kasier,'k','MarkerIndices',5:7:N,'LineWidth',1.3)
hold on
% Define the Fourier Transform using math calculation
w = 10;
f = -w:2*w/N:w-2*w/N;
S = 80*sinc(4*f);

% Truncation
S(1:threshold) = zeros(1,threshold);
S(N-threshold+1:N) = zeros(1,threshold);

% Gaussian reconstruction
for sigma = 1:1:150
    reconstruction_gaussian = G_filter(S,sigma);
    error_gau = abs(reconstruction_gaussian) - abs(s);
    points_Gau = find(abs(error_gau) < ep1);
    length_Gau(sigma,1) = length(points_Gau);
end
% Best sigma for Gaussian
Gaussian_best_sigma = median(find(length_Gau == max(length_Gau)));
H=zeros(1,N);
% Calculate using the transfer function
for i=1:N
    D=sqrt((i-1-floor(N/2))^2);
    H(1,i)=exp(-D^2/(2*Gaussian_best_sigma^2));
end
plot(xf,H,'-ko','MarkerIndices',1:7:N,'LineWidth',1.3)
hold on


% Define the Fourier Transform using math calculation
w = 10;
f = -w:2*w/N:w-2*w/N;
S = 80*sinc(4*f);

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
        error_ffre_exp = abs(reconstruction_exp) - abs(s);
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
    error_ffre_exp = abs(reconstruction_exp) - abs(s);
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
error_ffre_exp = abs(reconstruction_exp) - abs(s);
points_exp = find(abs(error_ffre_exp) < ep1);
% Save the performance of all framelen
result_exp(1,total_step+1-threshold) = length(points_exp);

plot(xf,exponentialfilter,'-k*','MarkerIndices',3:7:N,'LineWidth',1.3)
hold on


% cutf = 10-20*threshold/N;
Wn = 0.23*0.9;
order_save = zeros(1,50);
for order = 1:1:10
    [reconstructionButterW,~] = MyButter(order,Wn,S);
    error_butter = abs(reconstructionButterW) - abs(s);
    points_butter = find(abs(error_butter) < ep1);
    % Save the performance of all order
    order_save(order) = length(points_butter);
end
% Find best order
order_max = find(order_save == max(order_save),1 );
% Use best order to reconstruct
[reconstructionButterW,g_butter] =MyButter(order_max,Wn,S) ;
plot(xf,g_butter)







for Wn = 0.01:0.01:0.99
        order_save = zeros(1,50);
        for order = 1:1:10
            [reconstructionButterW,~] = MyButter(order,Wn,S);
            error_butter = abs(reconstructionButterW) - abs(s);
            points_butter = find(abs(error_butter) < ep1);
            % Save the performance of all order
            order_save(order) = length(points_butter);
        end
        % Find best order
        order_max = find(order_save == max(order_save),1 );
        % Use best order to reconstruct
        [reconstructionButterW,~] =MyButter(order_max,Wn,S) ;
        error_butter = abs(reconstructionButterW) - abs(s);
        points_butter = find(abs(error_butter) < ep1);
        % Save the performance of all 
        L0_ep_butter(round(100*Wn)) = length(points_butter);
        % Save the correspondent power
        order_cr(round(100*Wn)) = order_max;
    end
    re_butterW = find(L0_ep_butter == max(L0_ep_butter),1 );
    [reconstructionButterW,g_butter_grid] = MyButter(order_cr(re_butterW),re_butterW/100,S) ;
    error_butter = abs(reconstructionButterW) - abs(s);
    points_butter = find(abs(error_butter) < ep1);
    % Save the performance and the corresponding parameters
    result_butter(1,total_step+1-threshold) = length(points_butter);
    
    result_butter_order(1,total_step+1-threshold) = order_cr(re_butterW);
    result_butter_Wn(1,total_step+1-threshold) = re_butterW/100;
plot(xf,g_butter_grid)



xline(10-20*threshold/N,'-.b','LineWidth',1.3);
hold on
xline(-(10-20*threshold/N),'-.b','LineWidth',1.3);
legend('Hamming','Parzen window','Kasier window','Gaussian','exponential','Butterworth Wn = 0.2','Butterworth Grid search','Cut-Off Frequency')
xlabel('Frequency[Hz]')
ylabel('?')

