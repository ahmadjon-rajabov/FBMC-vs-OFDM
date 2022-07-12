%% OFDM
clear all;
close all;
clc

M = 16;                             %   16-QAM constellation
N = 2^20;                           %   number of symbols
no_of_ifft_points = 2^6;            %   points for the FFT/IFFT
no_of_fft_points = 2^6;             
block_size = no_of_fft_points;      %   number of symbols per ofdm channel
cp_len = ceil(0.1*block_size);      %   length of cyclic prefix
coordinates = [-3 -1 1 3];          
numSamplePerSymbol = 2;

tic
%% Transmitter side
%   Signal modulation
input = randsrc(1,N,coordinates) + 1i*randsrc(1,N,coordinates);
 
%   First step: obtain the number of OFDM channels (i.e. columns) that will exist after reshaping
num_ch = length(input)/block_size;
data_matrix = reshape(input, block_size, num_ch);

x = reshape(data_matrix, 1, (block_size * num_ch)); %input data in series form for eyediagram

%   Second: Create empty matix to put the IFFT'd data
cp_start = block_size-cp_len;
cp_end = block_size;
%   Third: Operate columnwise & do CP
for i=1:num_ch
    ifft_data_matrix(:,i) = ifft((data_matrix(:,i)),no_of_ifft_points);
    %   Compute and append Cyclic Prefix
    for j=1:cp_len
       actual_cp(j,i) = ifft_data_matrix(j+cp_start,i);
    end
    %   Append the CP to the existing block to create the actual OFDM block
    ifft_data(:,i) = vertcat(actual_cp(:,i),ifft_data_matrix(:,i));
end
%   Convert parallel to serial for transmission
[rows_ifft_data cols_ifft_data] = size(ifft_data);
len_ofdm_data = rows_ifft_data * cols_ifft_data;
%   Actual OFDM signal to be transmitted
ofdm_signal = reshape(ifft_data, 1, len_ofdm_data);


%%  Channel 
snr_dB = [0:20];  % here we use a loop to vary snr
errors = zeros(size(snr_dB));
for ii = 1:length(snr_dB)
    s = 1/sqrt(mean(abs(ofdm_signal).^2)); % 16-QAM normalization 
    n = 1/sqrt(2)*(randn(1,len_ofdm_data) + 1i*randn(1,len_ofdm_data)); %  normalized guassian noise
%   Pass the ofdm signal through the channel
    recvd_signal = s * ofdm_signal + 10^(-snr_dB(ii)/20) * n; % linear AWGN

%%  Receiver side
%   Convert Data back to "parallel" form to perform FFT
recvd_signal_matrix = reshape(recvd_signal,rows_ifft_data, cols_ifft_data);
%   Remove CP
recvd_signal_matrix(1:cp_len,:)=[];
%   Perform FFT
for i=1:cols_ifft_data
    fft_data_matrix(:,i) = fft(recvd_signal_matrix(:,i),no_of_fft_points);
end
%   Convert parallel to serial
y = reshape(fft_data_matrix, 1,(block_size*num_ch));
y = y./s;
%   Hard decision 
    y_re = real(y); % real part
    y_im = imag(y); % imaginary part
    out_re = y_re; out_im = y_im;
    out_re(y_re < -2) = -3;
    out_re(y_re > 2) =  3;
    out_re(y_re > -2 & y_re <= 0) = -1;
    out_re(y_re > 0 & y_re <= 2) =  1;
 
    out_im(y_im < -2) = -3;
    out_im(y_im > 2) =  3;
    out_im(y_im > -2 & y_im <= 0) = -1;
    out_im(y_im > 0 & y_im <= 2) =  1;
    out = out_re + 1i * out_im;

errors(ii) = length(find(input - out)); % calculate errors
end

%% Graph
scatterplot(input); title('Modulated data');
scatterplot(y); title('Received symbols before mapping');
scatterplot(out); title('Received symbols');

%to plot OFDM signal
plot(real(ofdm_signal)); xlabel('Time'); ylabel('Amplitude'); 
title('OFDM Signal');grid on;

%SER
Ser = errors / N;
formula_Ser = 3/2 * erfc(sqrt(0.1 * (10.^(snr_dB/10)))); %theoritical
figure
plot(snr_dB, log10(Ser),'-x'); hold on; %numerics
plot(snr_dB, log10(formula_Ser),'o'); %analytics
legend('numerics','analytics');
xlabel('SNR, dB');
ylabel('Lg(SER)');
title('Symbol Error Rate (a)');

%BER
Ber = Ser / 4 ; % a number of bits per symbol 16 QAM
formula_Ber = (1/4.0) * 3/2 * erfc(sqrt(0.1 * (10.^(snr_dB/10)))); %theoritical

figure 
plot(snr_dB, log10(Ber),'-x'); hold on; %numerics
plot(snr_dB, log10(formula_Ber),'o'); %analytics
legend('numerics','analytics');
xlabel('SNR, dB');
ylabel('Lg(BER)');
title('BER vs SNR (a)');

%Eye Diagram
eyediagram(x(1:1000), numSamplePerSymbol^2); %transmitter
eyediagram(y(1:1000), numSamplePerSymbol^2); %receiver

%Q-factor
Q = 20 * log10(sqrt(2) * (erfcinv(2 * (Ber))));
formula_Q = 20 * log10(sqrt(2) * (erfcinv(2 * (formula_Ber))));

figure 
plot(snr_dB, Q,'-x'); hold on; %numerics
plot(snr_dB, formula_Q,'o'); %analytics
title('Q factor vs SNR (a)');
ylabel('Q factor');
xlabel('SNR (dB)');
grid on
legend('numerics','analytics');

%Bandwidth efficiency
BW = no_of_fft_points/(no_of_fft_points + cp_len);
fprintf('Bandwidth efficiency is %.2f bit/s/Hz\n',BW)

toc
%%






