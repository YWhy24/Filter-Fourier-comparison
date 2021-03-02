clear all

width = 20;
N = width*20;
% Determining the x axis for all the signals 
x = -width/2:width/N:width/2-width/N;

% Define the signal
s = zeros(1, N);
t2 = (2*N/5)+1:(3*N/5);
s(t2) = 1;


% Define sigma range
number_of_sigma_tested = 150;
sigma_start = 1;
sigma_range = sigma_start:1:sigma_start+number_of_sigma_tested-1;


% The total number of truncations tested
total_step = 199;
% the truncation range. 
truncation_range = 1:1:total_step;
result = zeros(length(truncation_range),3);
w = 10;
f = -w:2*w/N:w-2*w/N;
S = 80*sinc(4*f);
% This is a waitbar so that the progress can be known
f2 = waitbar(0,'Processing');
% Special Fourier reconstruction
s_Fourier_special = fftshift(ifft(fftshift(S)));
% This is the standard for each truncation
ep1 =  median(abs(s-s_Fourier_special));
for threshold = truncation_range
%   update the waitbar
    waitbar(threshold/total_step)
% Fourier reconstruction
    S(1:threshold) = zeros(1,threshold);
    S(N-threshold+1:N) = zeros(1,threshold);
    s_Fourier_re = fftshift(ifft(fftshift(S)));
   
    L0_ep = zeros(number_of_sigma_tested,1);
    metrics = zeros(length(sigma_range),13);
    for sigma = sigma_range
        reconstruction = G_filter(S,sigma);
        error_s1 = reconstruction - s;
        ep0 = find(abs(error_s1) < ep1);
        metrics(sigma-sigma_start+1,:) =  metrics_calculation(s,reconstruction, s_Fourier_re);
        L0_ep(sigma) = length(ep0);
    end
%     metrics = metrics_normalization(s,S,metrics,s_Fourier_re,number_of_sigma_tested);
    p1 = median(find(L0_ep == max(L0_ep)));
    ppi = metrics(:,11);
    p2 = find(ppi == max(ppi));
%     beta = metrics(:,5);
%     p3 = find((beta(10:150,1) > 0.99),1);
    result(200-threshold,1) = p1;
    result(200-threshold,2) = p2;
%     result(threshold,3) = p3;
end
delete(f2)
%%
linew = 1.5;
linew_zi = 0.8;
FontS = 18;
plot(truncation_range/20, result(:,1),'b','LineWidth',linew )
hold on
plot(truncation_range/20, result(:,2),'r--','LineWidth',linew )
hold on
xlabel('Cut-Off Frequency [Hz]','FontSize',FontS ,'interpreter','latex')
ylabel('Best $\sigma$ ','FontSize',FontS ,'interpreter','latex')
legend('$\ell_\epsilon^0$' ,'PPI','FontSize',FontS,'interpreter','latex')


