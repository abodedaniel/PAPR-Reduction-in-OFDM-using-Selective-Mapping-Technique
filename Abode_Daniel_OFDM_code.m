%% Code on OFDM Exercise for Ph.D. Interview
%  Written by: Abode Daniel
%  Submitted:  Jeremy Nadal on 27th of Jan, 2021 
close all
clear all
clc

%% Declaring OFDM System Properties
M = 64;                        % 64-QAM constellation
k = log2(M);                   % number of bits per symbol
N = k*2^16;                    % number of binary sequence
num_sc = 256;                  % Number of Subcarrier
num_dsc = 16;                  % Number of data subcarrier
block_size = N/(num_dsc*k);    % number of symbols per subcarrier
cp_len = 64;                   % length of cyclic prefix

%% Transmitter side
ran = RandStream('swb2712');
input1 = randsrc(1,N,[0,1],ran);   % Binary Sequence
input = qammod(input1',M,'gray','InputType','bit','UnitAveragePower',true); %64QAM Modulation
% Allocating the 64QAM Symbol to data subcarrier
data_matrix = reshape(input, block_size, num_dsc);   
% ifft to generate 256 subcarriers
ifft_data_matrix = ifft(data_matrix',num_sc)';
% Compute and append Cyclic Prefix
cp_start = block_size-cp_len+1;
actual_cp = ifft_data_matrix(cp_start:end,:);
ifft_data = [actual_cp;ifft_data_matrix];
% Convert parallel to serial for transmission
[rows_ifft_data, cols_ifft_data]=size(ifft_data);
len_ofdm_data = rows_ifft_data*cols_ifft_data;
ofdm_signal = reshape(ifft_data, 1, len_ofdm_data); % Actual OFDM signal to be transmitted

%% Channel
EbNo = [-10:20]; % Eb/No 
errors=zeros(size(EbNo));
for ii = 1:length(EbNo)
snrdB = EbNo(ii) + 10*log10(k);  %Conver Eb/No to snr in dB
s = 1/sqrt(mean(abs(ofdm_signal).^2)); % Normalizer
n = 1/sqrt(2)*(randn(1,len_ofdm_data) + 1i*randn(1,len_ofdm_data)); % normalized guassian noise
% Pass the ofdm signal through the channel
recvd_signal = s*ofdm_signal+10^(-snrdB/20)*n; % linear AWGN

%% Receiver side
% Convert Data back to "parallel" form to perform FFT
recvd_signal_matrix = reshape(recvd_signal,rows_ifft_data, cols_ifft_data);
% Remove CP
recvd_signal_matrix(1:cp_len,:)=[];
% Perform FFT
fft_data_matrix = fft(recvd_signal_matrix',num_sc)';
fft_data_matrix(:,num_dsc+1:end) = [];  %Extract Data Symbols
% Convert parallel to serial and Normalize
y = reshape(fft_data_matrix, 1,[]);
y=y./s;

% Hard decision
out1 = qamdemod(y,M,'gray','OutputType','bit','UnitAveragePower',true);
out = reshape(out1,1,[]);

errors(ii) = length(find(input1- out)); % calculate errors
end
ber1 = (errors/length(input1));  %simulation ber
ber = berawgn(EbNo,'qam',M);     %theoretical ber of single carrier

%% BER plot
figure
semilogy(EbNo,ber1,'-x','linewidth',1.2);
hold on;
semilogy(EbNo,ber,'--','linewidth',1.2)
hold off
legend('simulation ber 16 of 256 subcarrier used','Theoretical ber of a single carrier','FontSize',15);
title('BER Performance','FontSize',15)
xlabel('Eb/No, dB','FontSize',15)
ylabel('Log10(BER)','FontSize',15)
ylim([0.5e-5 1])

%% Power Complementary Cumulative Distribution Function Derivation for OFDM
for i = 1:num_sc   %deriving the papr of all subcarriers
PAPR(i)= 10*log10(max(abs(ifft_data(:,i)).^2)/mean(abs(ifft_data(:,i))).^2);
end
[Y,X] = hist(PAPR,200);
figure
plot(X,1-cumsum(Y)/max(cumsum(Y)),'-b', 'LineWidth',1.2);
title('Power CCDF of OFDM','FontSize',15)
xlabel('Power above Average Power, dB','FontSize',15)
ylabel('Probability','FontSize',15)

%% Q4 Power Complementary Cumulative Distribution Function Derivation for Single Carrier
for i = 1:65536
PAPR2(i) = 10*log10(max(abs(input.^2))/mean(abs(input(i)).^2));
end
[Y1,X1] = hist(PAPR2,200);
figure
plot(X1,1-cumsum(Y1)/max(cumsum(Y1)),'-b', 'LineWidth',1.2);
title('Power CCDF of Single Carrier 64QAM','FontSize',15)
xlabel('Power above Average Power, dB','FontSize',15)
ylabel('Probability','FontSize',15)

%% Q4 Compare CCDF of OFDM and Single Carrier 64QAM
figure
plot(X1,1-cumsum(Y1)/max(cumsum(Y1)),'-b', 'LineWidth',1.2);
title('Power CCDF Comparison of Single Carrier 64QAM and OFDM','FontSize',15)
xlabel('Power above Average Power, dB','FontSize',15)
ylabel('Probability','FontSize',15)
hold on
plot(X,1-cumsum(Y)/max(cumsum(Y)), 'LineWidth',1.2);
legend('Power CCDF of 64QAM Single Carrier','Power CCDF of OFDM','FontSize',15)
hold off

%% Q5 Redesign of OFDM Transmitter using Selective Mapping Techniques
for j = 1:20                           %Performing the operation 40 times
B = randsrc(1,num_dsc,[1 -1 -1i 1i],ran);  %generate the phase sequence 
for i = 1:block_size
    d(i,:) = data_matrix(i,:).*B;      %Multiplying our 64QAM data matrix by the random Phase Sequence
end
sig = ifft(d',num_sc);                 %ifft
for i = 1:num_sc
    PAPR3(i)=10*log10(max(abs(sig(i,:)).^2)/mean(abs(sig(i,:))).^2); %calculating the PAPR
end
paprsum(j) = sum(PAPR3);                                             %summing the PAPR
% select signal with minimum PAPR
if(paprsum(j) == min(paprsum (1:j)))
    sigmin = sig;
end
end

%Calculating PAPR of selected signal
for i = 1:num_sc
    PAPR3(i)=10*log10(max(abs(sigmin(i,:)).^2)/mean(abs(sigmin(i,:))).^2);
end

figure
[Y3,X3] = hist(PAPR3,100);
plot(X3,1-cumsum(Y3)/max(cumsum(Y3)),'-.b', 'LineWidth',1.2);
title('Power CCDF Comparison of SLM Modified OFDM and OFDM','FontSize',15)
xlabel('Power above Average Power, dB','FontSize',15)
ylabel('Probability','FontSize',15)
xlim([9 12])
hold on
plot(X,1-cumsum(Y)/max(cumsum(Y)), 'LineWidth',1.2);
legend('Power CCDF of SLM Modified OFDM','Power CCDF of OFDM','FontSize',15)
hold off
