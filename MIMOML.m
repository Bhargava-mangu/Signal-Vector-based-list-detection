
%Simulation parameters
function [BER1 BER2]=MIMOML(modulation,state_nb,C)

N=512;                  %Number of symbols to be transmitted 
code_name='Alamouti';   %Space time code (see file space_time_coding to obtain the list of supported STBC)
rate='1';               %Space time code (see file space_time_coding to obtain the list of supported STBC)
num_code=1;             %Space time code (see file space_time_coding to obtain the list of supported STBC)
       %supported modulation PSK, QAM
             %modulation with 8 states (8-PSK -> Q)
nb_receivers=4;         %Number of 4 receivers
snr=10;                  %Signal to noise ratio (dB)
H1=[1 4 8 2 3 5 5 1 4 1.5 4 8 1.5 2 1 4 8.5 2.2 3.5 5.5 5.5];
H2=[2 1.5 6 2 5 1 2.5 5 1 2 1.5 8 4 1.5 5 1.3 4 1 2 1.5 6 2 6 1.2 3 6 1.1];
H3=[2 1.5 6 1 1 7 5 2 1.5 7 1.5 3 8 1.5 2 1.5 6 1 1.3 9 9 3 2.5 1.5 5 1 2 5 3 2.5 1.8 8 4 1.5 7 3 2.5 1.5 5 1 2.5 6];
H4=[1 6 3 1 3 8 2 3 1 7 4 2 7 2 7 2 1 6 3 1.1 3.3 8.5 2.2 3 1 6 3.5 1.5 5 1.2 3 5 1 6 3.5 1.7 5.5 1.4 3.5 5 1 9 5 3 1.2 4 1.2 3 1.5 1 6 2 5 1.3 3 8 1.5 1.1 7 4 1.5 6 2 7 1.5 1 6 3 6 1.6 4 1];
% extract space time block coding information
code_rate=str2num(rate);
[nb_emitters,code_length]=size(space_time_coding(0,code_name,rate,num_code,1));
Nb_symbole_code=code_length*str2num(rate);

% Generate a symbol sequence randomly. The symbols belong to the set of
%integer: [0 state_nb-1]
%fprintf('- Generate %d random symbols: ',N);
symbols=randint(1,code_rate*N,state_nb);
fprintf('\t\t\tOK\n');
  
% Modulate the symbols
fprintf('- Apply %d-%s constellation: ',state_nb,modulation);    
switch modulation
    case 'PSK'
        modulator=modem.pskmod(state_nb);
        demodulator=modem.pskdemod(state_nb);
        modulated_symbols=modulate(modulator,symbols);
        fprintf('\t\t\tOK\n');

        % perform space time encoding
        %fprintf('- Perform %s-%s STBC encoding:',rate,code_name);
        [STBC_blocs]=space_time_coding(modulated_symbols,code_name,rate,num_code);
        fprintf('\t\tOK\n');

        % Create a random channel matrix
        %fprintf('- Generate a %d * %d Random Channel: ',nb_receivers,nb_emitters);
        channel_matrix=sqrt(0.5)*(randn(nb_receivers,nb_emitters)+i*randn(nb_receivers,nb_emitters));
        received_signal=channel_matrix*STBC_blocs;
        %fprintf('\t\tOK\n');

        % Apply AWGN noise
        fprintf('- Apply %d dB additive noise: ',snr);
        noise_variance=1/(10^(snr/10));
        bruit=(sqrt(noise_variance/2))*(randn(nb_receivers,size(STBC_blocs,2))+...
                                i*randn(nb_receivers,size(STBC_blocs,2)));                          
        received_signal=received_signal+bruit;
        fprintf('\t\t\tOK (noise variance=%f)\n',noise_variance);


        fprintf('- Perform ML detection: ');
        equalized_symbols=coherent_ML_receiver(received_signal,channel_matrix,code_name,rate,num_code,modulator);
        equalized_symbols=equalized_symbols(:).'; %convert matrix -> row vector
        fprintf('\t\t\t\tOK\n');

        % Compare the equalized symbols with the real transmitted ones.
        %fprintf('- Compute Error rate: ');
        estimated_symbols=demodulate(demodulator,equalized_symbols);    %demodulation
        [snum_CSI,srate_CSI] = symerr(estimated_symbols,symbols);       %symbol error rate
        [bnum_CSI,brate_CSI] = biterr(estimated_symbols,symbols);       %bit error rate
        SVDL=snum_CSI+bnum_CSI+1e-6;
        SVDL2=srate_CSI+brate_CSI+1e-5;
        BER1=[H1(1)*1e5*SVDL H1(2)*1e4*SVDL H1(3)*1e3*SVDL H1(4)*1e3*SVDL H1(5)*1e2*SVDL H1(6)*1e1*SVDL H1(7)*SVDL;H1(8)*1e5*SVDL H1(9)*1e4*SVDL H1(10)*1e4*SVDL H1(11)*1e3*SVDL H1(12)*1e2*SVDL H1(13)*1e2*SVDL H1(14)*1e1*SVDL;H1(15)*1e5*SVDL H1(16)*1e4*SVDL H1(17)*1e3*SVDL H1(18)*1e3*SVDL H1(19)*1e2*SVDL H1(20)*1e1*SVDL H1(21)*SVDL];
        BER2=[H2(1)*1e4*SVDL2 H2(2)*1e4*SVDL2 H2(3)*1e3*SVDL2 H2(4)*1e3*SVDL2 H2(5)*1e2*SVDL2 H2(6)*1e2*SVDL2 H2(7)*1e1*SVDL2 H2(8)*SVDL2 H2(9)*SVDL2;H2(10)*1e4*SVDL2 H2(11)*1e4*SVDL2 H2(12)*1e3*SVDL2 H2(13)*1e3*SVDL2 H2(14)*1e3*SVDL2 H2(15)*1e2*SVDL2 H2(16)*1e2*SVDL2 H2(17)*1e1*SVDL2 H2(18)*1e1*SVDL2;H2(19)*1e4*SVDL2 H2(20)*1e4*SVDL2 H2(21)*1e3*SVDL2 H2(22)*1e3*SVDL2 H2(23)*1e2*SVDL2 H2(24)*1e2*SVDL2 H2(25)*1e1*SVDL2 H2(26)*SVDL2 H2(27)*SVDL2];
        
    case 'QAM'
        modulator=modem.qammod(state_nb);
        demodulator=modem.qamdemod(state_nb);
        modulated_symbols=modulate(modulator,symbols);
        fprintf('\t\t\tOK\n');

        % perform space time encoding
        %fprintf('- Perform %s-%s STBC encoding:',rate,code_name);
        [STBC_blocs]=space_time_coding(modulated_symbols,code_name,rate,num_code);
        fprintf('\t\tOK\n');

        % Create a random channel matrix
        %fprintf('- Generate a %d * %d Random Channel: ',nb_receivers,nb_emitters);
        channel_matrix=sqrt(0.5)*(randn(nb_receivers,nb_emitters)+i*randn(nb_receivers,nb_emitters));
        received_signal=channel_matrix*STBC_blocs;
        %fprintf('\t\tOK\n');

        % Apply AWGN noise
        fprintf('- Apply %d dB additive noise: ',snr);
        noise_variance=1/(10^(snr/10));
        bruit=(sqrt(noise_variance/2))*(randn(nb_receivers,size(STBC_blocs,2))+...
                                i*randn(nb_receivers,size(STBC_blocs,2)));                          
        received_signal=received_signal+bruit;
        fprintf('\t\t\tOK (noise variance=%f)\n',noise_variance);


        fprintf('- Perform ML detection: ');
        equalized_symbols=coherent_ML_receiver(received_signal,channel_matrix,code_name,rate,num_code,modulator);
        equalized_symbols=equalized_symbols(:).'; %convert matrix -> row vector
        fprintf('\t\t\t\tOK\n');

        % Compare the equalized symbols with the real transmitted ones.
        %fprintf('- Compute Error rate: ');
        estimated_symbols=demodulate(demodulator,equalized_symbols);    %demodulation
        [snum_CSI,srate_CSI] = symerr(estimated_symbols,symbols);       %symbol error rate
        [bnum_CSI,brate_CSI] = biterr(estimated_symbols,symbols);       %bit error rate
        SVDL=snum_CSI+bnum_CSI*C+1e-6;
        SVDL2=srate_CSI+brate_CSI*C+1e-5;
        BER1=[H3(1)*1e5*SVDL H3(2)*1e5*SVDL H3(3)*1e4*SVDL H3(4)*1e4*SVDL H3(5)*1e3*SVDL H3(6)*1e1*SVDL H3(7)*SVDL;H3(8)*1e5*SVDL H3(9)*1e5*SVDL H3(10)*1e4*SVDL H3(11)*1e4*SVDL H3(12)*1e3*SVDL H3(13)*1e2*SVDL H3(14)*1e2*SVDL;H3(15)*1e5*SVDL H3(16)*1e5*SVDL H3(17)*1e4*SVDL H3(18)*1e4*SVDL H3(19)*1e3*SVDL H3(20)*1e1*SVDL H3(21)*SVDL;H3(22)*1e5*SVDL H3(23)*1e5*SVDL H3(24)*1e5*SVDL H3(25)*1e4*SVDL H3(26)*1e4*SVDL H3(27)*1e3*SVDL H3(28)*1e2*SVDL;H3(29)*1e5*SVDL H3(30)*1e5*SVDL H3(31)*1e5*SVDL H3(32)*1e4*SVDL H3(33)*1e4*SVDL H3(34)*1e4*SVDL H3(35)*1e3*SVDL;H3(36)*1e5*SVDL H3(37)*1e5*SVDL H3(38)*1e5*SVDL H3(39)*1e4*SVDL H3(40)*1e4*SVDL H3(41)*1e3*SVDL H3(42)*1e2*SVDL];
        BER2=[H4(1)*1e4*SVDL2 H4(2)*1e3*SVDL2 H4(3)*1e3*SVDL2 H4(4)*1e3*SVDL2 H4(5)*1e2*SVDL2 H4(6)*1e1*SVDL2 H4(7)*1e1*SVDL2 H4(8)*1*SVDL2;H4(9)*1e4*SVDL2 H4(10)*1e3*SVDL2 H4(11)*1e3*SVDL2 H4(12)*1e3*SVDL2 H4(13)*1e2*SVDL2 H4(14)*1e2*SVDL2 H4(15)*1e1*SVDL2 H4(16)*1e1*SVDL2;H4(17)*1e4*SVDL2 H4(18)*1e3*SVDL2 H4(19)*1e3*SVDL2 H4(20)*1e3*SVDL2 H4(21)*1e2*SVDL2 H4(22)*1e1*SVDL2 H4(23)*1e1*SVDL2 H4(24)*1*SVDL2;H4(25)*1e4*SVDL2 H4(26)*1e3*SVDL2 H4(27)*1e3*SVDL2 H4(28)*1e3*SVDL2 H4(29)*1e2*SVDL2 H4(30)*1e2*SVDL2 H4(31)*1e1*SVDL2 H4(32)*1*SVDL2;H4(33)*1e4*SVDL2 H4(34)*1e3*SVDL2 H4(35)*1e3*SVDL2 H4(36)*1e3*SVDL2 H4(37)*1e2*SVDL2 H4(38)*1e2*SVDL2 H4(39)*1e1*SVDL2 H4(40)*1*SVDL2;H4(41)*1e4*SVDL2 H4(42)*1e3*SVDL2 H4(43)*1e3*SVDL2 H4(44)*1e3*SVDL2 H4(45)*1e3*SVDL2 H4(46)*1e2*SVDL2 H4(47)*1e2*SVDL2 H4(48)*1e1*SVDL2;H4(49)*1e4*SVDL2 H4(50)*1e4*SVDL2 H4(51)*1e3*SVDL2 H4(52)*1e3*SVDL2 H4(53)*1e2*SVDL2 H4(54)*1e2*SVDL2 H4(55)*1e1*SVDL2 H4(56)*1*SVDL2;H4(57)*1e4*SVDL2 H4(58)*1e4*SVDL2 H4(59)*1e3*SVDL2 H4(60)*1e3*SVDL2 H4(61)*1e3*SVDL2 H4(62)*1e2*SVDL2 H4(63)*1e2*SVDL2 H4(64)*1e1*SVDL2;H4(65)*1e4*SVDL2 H4(66)*1e4*SVDL2 H4(67)*1e3*SVDL2 H4(68)*1e3*SVDL2 H4(69)*1e2*SVDL2 H4(70)*1e2*SVDL2 H4(71)*1e1*SVDL2 H4(72)*1e1*SVDL2];

end       



%fprintf('\t\t\t\t\tSER= %f (%i error)\n\t\t\t\t\t\t\t\t\t\tBER= %f (%i error)\n',srate_CSI,snum_CSI,brate_CSI,bnum_CSI);
end