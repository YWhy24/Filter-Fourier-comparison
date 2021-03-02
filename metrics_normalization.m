function metrics_f = metrics_normalization(s,S,metrics,S_Fourier,number_of_sigma_tested)

N = length(s);
metrics_f = metrics;
                
min_mse = immse(s,S_Fourier );
max_mse = immse(s,G_filter(S,0.00001));

min_var_1 = var(  S_Fourier  );
max_var_1 = var(  G_filter(S,0.00001) );
    
min_entropy_s1 = pentropy(abs( G_filter(S,0.00001) ),N,'Instantaneous',false); % Calculate the entropy
max_entropy_s1 = pentropy(abs( S_Fourier ),N,'Instantaneous',false); % Calculate the entropy
   
min_psnr_s1 = 10*log10(max(abs( G_filter(S,0.00001)))^2 / max_mse ); 
max_psnrs1 = 10*log10(max(abs( S_Fourier))^2 / min_mse ); 

min_beta_s1 = 0;
max_beta_s1 = beta_calculation(s,S_Fourier );

min_var_2 = var( abs(s -  S_Fourier ) );
max_var_2 = var( abs( s - G_filter(S,0.00001)) );

min_max_error = max(abs( S_Fourier-s));
max_max_error = max(abs( s - G_filter(S,0.00001)));

min_snr = 10*log10(sum(s.^2)/sum((abs(S_Fourier-s)).^2));
max_snr = 10*log10(sum(s.^2)/sum((abs(G_filter(S,0.00001)-s)).^2));

max_energy = sum(S_Fourier.^2);
min_energy = sum(G_filter(S,0.00001).^2);
% 
min_hfen = hfen_calc(S_Fourier,s);
max_hfen =  hfen_calc(G_filter(S,0.00001),s);

min_psr = mean(sqrt(fft(G_filter(S,0.00001)).^2 ./(fft(S_Fourier)).^2));
max_psr = 1;


for a = 1:1:number_of_sigma_tested
    metrics_f(a,1) = (metrics_f(a,1)- min_mse)/(max_mse-min_mse);
end

for a = 1:1:number_of_sigma_tested
    metrics_f(a,2) = (metrics_f(a,2)- min_var_1)/(max_var_1-min_var_1);
end

for a = 1:number_of_sigma_tested
    metrics_f(a,3) = (metrics(a,3)-min_psnr_s1)  /  (max_psnrs1-min_psnr_s1);
end

for a = 1:number_of_sigma_tested
    metrics_f(a,4) = (metrics(a,4)-min_entropy_s1)  /  (max_entropy_s1-min_entropy_s1);
end

for a = 1:number_of_sigma_tested
    metrics_f(a,5) = (metrics(a,5)-min_beta_s1)  /  (max_beta_s1-min_beta_s1);
end

for a = 1:number_of_sigma_tested
    metrics_f(a,6) = (metrics(a,6)-min_var_2)  /  (max_var_2-min_var_2);
end

for a = 1:number_of_sigma_tested
    metrics_f(a,7) = (metrics(a,7)-min_max_error)  /  (max_max_error-min_max_error);
end

for a = 1:number_of_sigma_tested
    metrics_f(a,8) = (metrics(a,8)-min_snr)  /  (max_snr-min_snr);
end

for a = 1:number_of_sigma_tested
    metrics_f(a,9) = (metrics(a,9)-min_energy)  /  (max_energy-min_energy);
end

for a = 1:number_of_sigma_tested
    metrics_f(a,12) = (metrics(a,12)-min_hfen)  /  (max_hfen-min_hfen);
end

for a = 1:number_of_sigma_tested
    metrics_f(a,13) = (metrics(a,13)-min_psr)  /  (max_psr-min_psr);
end

end

