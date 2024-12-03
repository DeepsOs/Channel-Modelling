clear ;
close all;
clc;

N = 1:3; % No. of links
M = 1000000; % No. of samples

Eb_N0_dB = 0:5:40; %average SNR
   Eb_N0 = 10.^(Eb_N0_dB./10); %convert SNR to power

Eb_N0_th_dB = 5; % threshold SNR
   Eb_N0_th = 10.^(Eb_N0_th_dB./10); % convert SNR threshold to dB

m = 1;  % fading severity parameter

%---------------------Analytical------------------------

C = (m ./ Eb_N0);  % constant C=m/avg SNR
z = C .* Eb_N0_th;
g = igamma(m, z)/gamma(m);
gg= 1 - g;

for ii = 1:length(N)   
Pout_RF(ii,:) = (gg).^N(1,ii) ; % outage probability for MISO RF (Analytical)
end


%--------------------Simulations-------------------------
               
  %-------------N = 1 link------------------------------------------  
  
 for n = 1:3
 h(n,:) = sqrt(gamrnd(m, 1/m, 1, M) ); % Nakagami-m fading channel
 end
 
 
for k = 1: length(Eb_N0)
    
   instt_SNR(k,:) = Eb_N0(1,k).*((abs(h(1,:))).^2) ; % insttSNR = avgSNR * |h|^2
    
       
 for  jj = 1: size(instt_SNR,2)
   if  instt_SNR(k,jj) <  Eb_N0_th % checking Pout(instt SNR < threshold SNR)
      count(k,jj) = instt_SNR(k,jj);
       else
      count(k,jj) = 0; 
      end
 end
        nErr(:,k) = size(find(count(k,:)),2);  % counting the errrors
end
       P_out(1,:) =  nErr./M;  % outage probability for MISO N=1 (simulated)
   
 
    
  %------------------- N = 2 links --------------------------------      
  
  for kk = 1: length(Eb_N0)
             
     instt_SNR_1(kk,:) = Eb_N0(1,kk).*((abs(h(1,:))).^2) ; % insttSNR = avgSNR * |h|^2
     instt_SNR_2(kk,:) = Eb_N0(1,kk).*((abs(h(2,:))).^2) ;
     
  for s = 1:size(instt_SNR_1,2)
    
    if ((instt_SNR_1(kk,s)) < (instt_SNR_2(kk,s))) % selecting link with maximum instt SNR
       
    max_instt_SNR(kk,s) = instt_SNR_2(kk,s);
        
    else
    max_instt_SNR(kk,s) = instt_SNR_1(kk,s);
    end
 
  end
   
 for jjj = 1: size(max_instt_SNR,2)
     
   if  max_instt_SNR(kk,jjj) <  Eb_N0_th % checking Pout(instt SNR < threshold SNR)
       count(kk,jjj) = max_instt_SNR(kk,jjj);
       else
       count(kk,jjj) = 0; 
       end
 end
          nErr(:,kk) = size(find(count(kk,:)),2);

  end
          P_out(2,:) =  nErr./M;  % outage probability for MISO N=2 (simulated)
        
       
%----------------------- N = 3 links ---------------------------

 for kkk = 1: length(Eb_N0)
            
   instt_SNR_11(kkk,:) = Eb_N0(1,kkk).*((abs(h(1,:))).^2) ; % insttSNR = avgSNR * |h|^2
   
   instt_SNR_22(kkk,:) = Eb_N0(1,kkk).*((abs(h(2,:))).^2) ;    
      
   instt_SNR_33(kkk,:) = Eb_N0(1,kkk).*((abs(h(3,:))).^2) ; 
       
      
  for ss = 1:size(instt_SNR_11,2)
     if ((instt_SNR_11(kkk,ss)) > (instt_SNR_22(kkk,ss)) &&...
        ((instt_SNR_11(kkk,ss)) > (instt_SNR_33(kkk,ss)))) % selecting link with maximum instt SNR
            
   max_instt_SNR(kkk,ss) = instt_SNR_11(kkk,ss);
           
    elseif((instt_SNR_22(kkk,ss)) > (instt_SNR_11(kkk,ss)) &&...
          ((instt_SNR_22(kkk,ss)) > (instt_SNR_33(kkk,ss))))
   max_instt_SNR(kkk,ss) = instt_SNR_22(kkk,ss);
           
     else
   max_instt_SNR(kkk,ss) = instt_SNR_33(kkk,ss);
       end
       end
       
  for  jjjj = 1: size(max_instt_SNR,2)
     
  if(  max_instt_SNR(kkk,jjjj) <  Eb_N0_th )% checking Pout(instt SNR < threshold SNR)
    count(kkk,jjjj) = max_instt_SNR(kkk,jjjj);
       else
    count(kkk,jjjj) = 0; 
       end
  end
 
        nErr(:,kkk) = size(find(count(kkk,:)),2);

 end
         P_out(3,:) =  nErr./M ;  % outage probability for MISO N=3 (simulated)
  

figure

semilogy(Eb_N0_dB,Pout_RF,'Linewidth',2);
grid on
hold on
semilogy(Eb_N0_dB, P_out,'x','Linewidth',2);
hold on
xlabel('RF average SNR (Es/No) / dB');
ylabel(' Outage probability');
legend('RF link N = 1 - Analytical','RF link N = 2- Analytical','RF link N = 3 - Analytical','RF link N = 1 - Simulation', 'RF link N = 2 - Simulation','RF link N = 3 - Simulation','Location','SouthWest');
hold on

