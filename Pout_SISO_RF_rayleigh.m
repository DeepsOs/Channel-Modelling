clear;
close all;
clc;

N = 1; % no. of links
M = 100000; % no. of samples

Eb_N0_dB = 0:5:40 ; % multiple Eb/N0 values
   Eb_N0 = 10.^(Eb_N0_dB./10);   % convert dB to power

Eb_N0_th_dB = [5  10 15]; % SNR Threshold
   Eb_N0_th = 10.^(Eb_N0_th_dB./10); % convert dB to power
   
m = 1; % fading severity parameter

%--------------------------Analytical------------------------------

for i = 1:length(Eb_N0_th_dB)
    C = m ./ Eb_N0;  
    z = C .* Eb_N0_th(1,i); 
  % g = igamma(m, z)/gamma(m); % incomplete gamma function
 % gg = 1 - g;
gammainc = 1 - (igamma(m, z)/gamma(m));  % lower incomplete gamma function
Pout_RF(i,:) = (gammainc).^N ; % outage probability for SISO (analytical)
end


%-------------Simulations------------           

h = 1/sqrt(2)*[randn(1,M) + j*randn(1,M)]; % rayleigh fading channel
 
 for k = 1: length(Eb_N0)
     
        instt_SNR(k,:) = Eb_N0(k).*((abs(h)).^2) ; % insttSNR = avgSNR * |h|^2
     
     
 for i = 1:length(Eb_N0_th_dB)
      
 for  j = 1: size(instt_SNR,2)
       
       if  instt_SNR(k,j) < Eb_N0_th(1,i) % checking Pout(instt SNR < threshold SNR)
           count(i,j) = 1;
       else
           count(i,j) = 0; 
       end
 end
 
            nErr(:,i) = size(find(count(i,:)),2);
 end

        P_out(:,k)=  nErr/M;  % outage probability for SISO (simulated)
    end 



figure
semilogy(Eb_N0_dB,Pout_RF,'Linewidth',2);
grid on
hold on
semilogy(Eb_N0_dB, P_out,'x','Linewidth',2);
hold on
xlabel('RF average SNR (Es/No) / dB');
ylabel(' Outage probability');
legend('RF link avg_ snr = 5 dB- Analytical','RF link avg_ snr = 10 dB- Analytical','RF link avg_ snr = 15 dB- Analytical','RF link avg_ snr = 5 dB- Simulation','RF link avg_ snr = 10 dB- Simulation','RF link avg_ snr = 15 dB- Simulation','Location','Southwest');
hold on
