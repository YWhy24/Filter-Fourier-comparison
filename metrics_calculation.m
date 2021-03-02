function [metrics] = metrics_calculation(signal,reconstruction,Fourier_reconstruction)

    N = length(signal);
    metrics = zeros(1,11);
    
    mse = immse( signal, reconstruction ); % Calculate the mean squared error of the first signal
    variance = var( reconstruction ); % Calculate the variance of the first signal
    psnr = 10*log10(max(abs( reconstruction))^2 / mse); % Calculate the PSNR   
    entropy = pentropy(abs( reconstruction ),N,'Instantaneous',false); % Calculate the entropy
    beta = beta_calculation(signal,reconstruction); % Calaulate beta
    e = abs(reconstruction - signal);
    variance_2 = var(e );
    max_error = max(abs(reconstruction-signal));
    SNR = 10*log10(sum(abs(signal).^2)/sum((abs(reconstruction-signal)).^2));
    Energy = sum(abs(reconstruction).^2);
    correlation = correlationcalc(signal,reconstruction);
    % ?
    ppi = 1/(mean(sqrt(  abs(reconstruction).^2  ./  abs(Fourier_reconstruction).^2   )));
    hfen = hfen_calc(reconstruction,signal);
    % psr
    Actual_psr = mean(sqrt(fft(reconstruction).^2 ./(fft(Fourier_reconstruction)).^2));
    metrics(1) = mse;
    metrics(2) = variance;
    metrics(3) = psnr;
    metrics(4) = entropy;
    metrics(5) = beta;
    metrics(6) = variance_2; 
    metrics(7) = max_error;
    metrics(8) = SNR;
    metrics(9) = Energy;
    metrics(10) = correlation;
    metrics(11) = ppi;
    metrics(12) = hfen;
    metrics(13) = Actual_psr ;
end

