% This document is written by Yue Wang 
% Finalized on the 15th of August, 2020
% This is a filtered Fourier reconstruction function 
% Using the Gaussian function

function reconstruction = G_filter(S,sigma)
%The length of the filter
N = length(S);
%Initialization
H=zeros(1,N);
% Calculate using the transfer function
for i=1:N
    D=sqrt((i-1-floor(N/2))^2);
    H(1,i)=exp(-D^2/(2*sigma^2));
end
%Create the filtered reconstruction
reconstruction = fftshift(ifft(fftshift(S.*H)));

end

