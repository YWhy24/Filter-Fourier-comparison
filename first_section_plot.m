clear all
width = 20;
N = width*20;
% Determining the x axis for all the signals 
x = -width/2:width/N:width/2-width/N;
% Define the signal
s = zeros(1, N);
t2 = (2*N/5)+1:(3*N/5);
s(t2) = 1;
threshold = 60;
w = 10;
f = -w:2*w/N:w-2*w/N;
S = 80*sinc(4*f);



S_for_comparison = fftshift(fft(fftshift(s)));

linew = 0.8;
linew_zi = 0.8;
FontS = 18;

% 
% figure(1)% The original signal's fft
% 
% plot(f,1/20*abs(S_for_comparison),'--r','LineWidth',linew )
% hold on
% plot(f,1/20*abs(S),'b','LineWidth',linew )
% xlabel('Frequency f(Hz)','FontSize',FontS )
% ylabel('Value |G(f)|','FontSize',FontS )
% ylim([0,4.4])
% legend('Calculated using signal samples','Calculated using mathematical function','Location','northwest','FontSize',FontS )
% axes('position',[.6 .40 .25 .25])
% box on % put box around new pair of axes
% indexOfInterest = (f<10) & (f > 8); % range of t near perturbation
% plot(f(indexOfInterest),abs(S_for_comparison(indexOfInterest)),'--r','Linewidth',linew_zi) % plot on new axes
% hold on
% plot(f(indexOfInterest),abs(S(indexOfInterest)),'b','LineWidth',linew_zi )
% axis tight
% ylim([0,1.1])
% title('Zoomed In','FontSize',FontS )
% hold off

% Fourier reconstruction
S(1:threshold) = zeros(1,threshold);
S(N-threshold+1:N) = zeros(1,threshold);
s_Fourier_re = fftshift(ifft(fftshift(S)));


figure(2)% The original signal
plot(x,s,'--r','Linewidth',linew)

ylim([-0.1,1.1])
set(gca,'FontSize',FontS)
hold on 
plot(x,real(s_Fourier_re),'b','Linewidth',linew)
xlabel('Time(t)','FontSize',FontS,'interpreter','latex' )
ylabel('g(t)','FontSize',FontS ,'interpreter','latex')
legend('Example Signal','Fourier Reconstruction','Location','northwest','FontSize',FontS-3,'interpreter','latex' )
axes('position',[.65 .50 .2 .2])
box on % put box around new pair of axes
indexOfInterest = (x<2) & (x>0); % range of t near perturbation
plot(x(indexOfInterest),abs(s(indexOfInterest)),'--r','Linewidth',linew_zi) % plot on new axes
set(gca,'FontSize',FontS-3)
hold on
plot(x(indexOfInterest),abs(s_Fourier_re(indexOfInterest)),'b','Linewidth',linew_zi)
axis tight
xlim([1.2, 2.2])
ylim([0.8,1.1])
title('Zoomed In','FontSize',20,'interpreter','latex' )
hold off

% sigma definition
number_of_sigma_tested = 150;
sigma_start = 1;
sigma_range = sigma_start:1:sigma_start+number_of_sigma_tested-1;

ep1 =  abs(median(s-s_Fourier_re));
L0_ep = zeros(number_of_sigma_tested,1);
metrics = zeros(length(sigma_range),13);
for sigma = sigma_range
    reconstruction = G_filter(S,sigma);
    error_s1 = abs(reconstruction) - abs(s);
    p1 = find(abs(error_s1) < ep1);
    metrics(sigma-sigma_start+1,:) =  metrics_calculation(s,reconstruction,s_Fourier_re );
    L0_ep(sigma) = length(p1);
end

metrics = metrics_normalization(s,S,metrics,s_Fourier_re,number_of_sigma_tested);
re1 = find(L0_ep == max(L0_ep) );
% figure(3)
linew1 = 1.5;
% plot(sigma_range,real(metrics(:,11)),'Linewidth',1.2)
% 
% figure(3)
% plot(sigma_range,metrics(:,1:2),'Linewidth',1.2)
% hold on
% plot(sigma_range,metrics(:,3:4),'Linewidth',1.2)
% hold on
% plot(sigma_range,real(metrics(:,6:8)),'Linewidth',1.2)
% hold on
% plot(sigma_range,real(metrics(:,9:10)),'Linewidth',1.2)
% hold on
% plot(sigma_range,real(metrics(:,12)),'Linewidth',1.2)
% hold on
% plot(sigma_range,real(metrics(:,13)),'Linewidth',1.2)
% set(gca,'FontSize',FontS)
% hold on
% % xlim([1.2, 2.2])
% ylim([-0.1,1.1])
% xlabel('The \sigma Parameter of Gaussian Filter','FontSize',FontS )
% ylabel('Metric','FontSize',FontS )
% legend('MSE','Variance','PSNR','Entropy','Error Variance','Max Error','SNR','Energy','Correlation','HFEN','PSR','Location','east','FontSize',FontS-0.2 )
% % plot(sigma_range,(metrics(:,11)),'Linewidth',1.2)

% 
% figure(4)
% comb = (20*metrics(:,1)+metrics(:,3))/21;
% plot(sigma_range,comb,'b','Linewidth',1.1)
% set(gca,'FontSize',FontS)
% xlabel('The \sigma Parameter of Gaussian Filter','FontSize',FontS-4 )
% ylabel('Normalized Combination Metric Value','FontSize',FontS-4 )
% hold on
% % p1 = find(comb == min(comb));
% % xline(min(p1),'b','The \sigma with lowest possible value','Linewidth',1.1)
% axes('position',[.6 .40 .25 .25])
% box on % put box around new pair of axes
% indexOfInterest = (sigma_range<90) & (sigma_range > 70); % range of t near perturbation
% plot(sigma_range(indexOfInterest),comb(indexOfInterest),'b','Linewidth',linew_zi) % plot on new axes
% set(gca,'FontSize',FontS)
% title('Zoomed In','FontSize',FontS )
% axis tight
% ylim([0.04729,0.04740])
% % 
figure(5)
beta = metrics(:,5);
plot(sigma_range,beta,'b','Linewidth',2)
set(gca,'FontSize',FontS)
xlabel('The \sigma Parameter of Gaussian Filter','FontSize',FontS )
ylabel('Normalized \beta Metric Value','FontSize',FontS )
hold on
% p1 = find(beta == max(beta((30:150),:)));
% xline(min(p1),'b','The \sigma with lowest possible value','Linewidth',1.1)
axes('position',[.6 .40 .25 .25])
box on % put box around new pair of axes
indexOfInterest = (sigma_range<80) & (sigma_range > 50); % range of t near perturbation
plot(sigma_range(indexOfInterest),beta(indexOfInterest),'b','Linewidth',2) % plot on new axes
set(gca,'FontSize',FontS)
title('Zoomed In','FontSize',FontS )
axis tight
ylim([0.998,1.006])
% % 
% % 
figure(6)


yyaxis left
plot(sigma_range,L0_ep,'b','Linewidth',1.1)
set(gca,'FontSize',FontS)
xlabel('The \sigma Parameter of Gaussian Filter','FontSize',FontS )
ylabel('$\ell_\epsilon^0$','FontSize',FontS ,'interpreter','latex')
hold on
p1 = find(L0_ep == max(L0_ep));
% xline(min(p1),'b','The \sigma with lowest possible value','Linewidth',1.1)
% axes('position',[.4 .28 .25 .25])
% box on % put box around new pair of axes
% indexOfInterest = (sigma_range<65) & (sigma_range > 35); % range of t near perturbation
% plot(sigma_range(indexOfInterest),L0_ep(indexOfInterest),'b','Linewidth',linew_zi) % plot on new axes
% set(gca,'FontSize',FontS)
% title('Zoomed In','FontSize',FontS )
axis tight
yyaxis right
plot(sigma_range,real(metrics(:,11)),'r--','Linewidth',1.5)
ylabel('PPI','FontSize',FontS)
legend('$\ell_\epsilon^0$','PPI','FontSize',FontS ,'interpreter','latex')