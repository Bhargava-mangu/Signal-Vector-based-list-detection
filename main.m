clc;clear all;close all

fprintf('\n+-----------------------------------+\n');
fprintf('|   Simulate a MIMO communication   |\n');
fprintf('+-----------------------------------+\n\n');

nt_V = [4 4 4 5 8];
nr_V = [3 3 4 5 6];

N0 = 1e-4;
B  = 1;
Iteration = 1e4; % must be grater than 1e2

SNR_V_db = [2:4:26];
SNR_V    = 10.^(SNR_V_db/10);
modulation1='PSK';
state_nb1=8;
modulation2='QAM';
state_nb2=16;

for(k = 1 : 5)
    fprintf('|   Apply a SVD  |\n');

    nt = nt_V(k);
    nr = nr_V(k);
    for(i = 1 : length(SNR_V))
        Pt = N0 * SNR_V(i);
        for(j = 1 : Iteration)
            H = random('rayleigh',1,nr,nt);
            [S V D] = svd(H);
            landas(:,j)  = diag(V);
            [Capacity(i,j) PowerAllo] = WaterFilling_alg(Pt,landas(:,j),B,N0);
        end
    end
    C=sum(sum(Capacity));
    [BER1 BER2]=MIMOML(modulation1,state_nb1,C);
    [BER3 BER4]=MIMOML(modulation2,state_nb2,C);
    %BERSVD
    clear landas
end

snr1=2:2:14;
snr3=2:4:26;
snr2=2:2:18;
snr4=2:2:16;

figure
semilogy(snr2,BER2(1,:),'-ok',snr2,BER2(2,:),'-ks',snr2,BER2(3,:),'-kd',snr1,BER1(1,:),'--ok',snr1,BER1(2,:),'--ks',snr1,BER1(3,:),'--kd','LineWidth',2);
xlabel('SNR(dB)');ylabel('BER')
legend('ML, Nr=4','SVD, Nr=4','List SVD, L=2, Nr=4','ML, Nr=6','SVD, Nr=6','List SVD, L=2, Nr=6','Location','SouthWest');
title('BER performance curves with SVD, list SVD, and ML detection when Nt = 8 and X is 8-PSK');

figure
semilogy(snr3,BER3(1,:),'-ok',snr3,BER3(2,:),'-ks',snr3,BER3(3,:),'-kd',snr3,BER3(4,:),'--ok',snr3,BER3(5,:),'--ks',snr3,BER3(6,:),'--kd','LineWidth',2);
xlabel('SNR(dB)');ylabel('BER')
legend('ML, Nr=3','SVD, Nr=3','List SVD, L=2, Nr=3','ML, Nr=2','SVD, Nr=2','List SVD, L=2, Nr=2','Location','SouthWest');
title('BER performance curves with SVD, list SVD, and ML detection when Nt = 4 and X is 16-QAM');
axis([2 26 1e-6 1])

figure
semilogy(snr4,BER4(1,:),'-ok',snr4,BER4(2,:),'-ks',snr4,BER4(3,:),'-kd',snr4,BER4(4,:),'--ok',snr4,BER4(5,:),'--ks',snr4,BER4(6,:),'--kd','LineWidth',2);
xlabel('SNR(dB)');ylabel('BER')
legend('ML, p=0.2','SVD, p=0.2','List SVD, L=2, p=0.2','ML, p=0.6','SVD, p=0.6','List SVD, L=2, p=0.6','Location','SouthWest');
title('BER performance curves with SVD, list SVD, and ML detection when Nt = Nr = 4 and X is 8-PSK, in the correlated fading case');
axis([2 16 1e-6 1])

figure
semilogy(snr4,BER4(7,:),'-ok',snr4,BER4(8,:),'-ks',snr4,BER4(9,:),'-kd','LineWidth',2);
xlabel('SNR(dB)');ylabel('BER')
legend('ML','SVD','SVD, L=2','Location','SouthWest');
title('BER performance curves with SVD, list SVD, and ML detection when Nt = 5,Nr = 4 and X is 8-PSK, in the GSM');
axis([2 16 1e-6 1])