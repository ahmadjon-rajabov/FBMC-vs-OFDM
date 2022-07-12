clear all
clc
close 
 
s = rng(211);   
 
numFFT = 1024;           % Number of FFT points
numGuards = 212;         % Guard bands on both sides
K = 4;                   % Overlapping symbols, one of 2, 3, or 4
numSymbols = 2^10;       % Simulation length in symbols
bitsPerSubCarrier = 4;   % 2: 4QAM, 4: 16QAM, 6: 64QAM, 8: 256QAM
snrdB = [0:20];          % SNR in dB
num_samples_each_trace = 4; %number of samples for each trace for eye
    %diagram

tic 
% Prototype filter
switch K
    case 2
        HkOneSided = sqrt(2)/2;
    case 3
        HkOneSided = [0.911438 0.411438];
    case 4
        HkOneSided = [0.971960 sqrt(2)/2 0.235147];
    otherwise
        return
end
% Build symmetric filter
Hk = [fliplr(HkOneSided) 1 HkOneSided];
 
% Transmit-end processing
% Initialize arrays
L = numFFT-2*numGuards;  	% Number of complex symbols per OFDM symbol
KF = K*numFFT;
KL = K*L;
dataSubCar = zeros(L, 1);
dataSubCarUp = zeros(KL, 1);
 
sumFBMCSpec = zeros(KF*2, 1);
sumOFDMSpec = zeros(numFFT*2, 1);
 
numBits = bitsPerSubCarrier*L/2;    % account for oversampling by 2
inpData = zeros(numBits, numSymbols);
rxBits = zeros(numBits, numSymbols);
txSigAll = complex(zeros(KF, numSymbols));
symBuf = complex(zeros(2*KF, 1));
 
% Loop over symbols
for symIdx = 1:numSymbols
 
%    Generate mapped symbol data
      inpData(:, symIdx) = randi([0 1], numBits, 1);
      modData = qammod(inpData(:, symIdx), 2^bitsPerSubCarrier, ...
         'InputType', 'Bit', 'UnitAveragePower', true);
 
    % OQAM Modulator: alternate real and imaginary parts
    %The real and imaginary parts of a complex data symbol 
    %are not transmitted simultaneously, as the imaginary part is delayed
    %by half the symbol duration.
    if rem(symIdx,2)==1     % Odd symbols
        dataSubCar(1:2:L) = real(modData);
        dataSubCar(2:2:L) = 1i*imag(modData);
    else                    % Even symbols
        dataSubCar(1:2:L) = 1i*imag(modData);
        dataSubCar(2:2:L) = real(modData);
    end
 
    % Upsample by K, pad with guards, and filter with the prototype filter
    dataSubCarUp(1:K:end) = dataSubCar;
    dataBitsUpPad = [zeros(numGuards*K,1); dataSubCarUp; zeros(numGuards*K,1)];
    X1 = filter(Hk, 1, dataBitsUpPad);
    % Remove 1/2 filter length delay
    X = [X1(K:end); zeros(K-1,1)];
 
    % Compute IFFT of length KF for the transmitted symbol
    txSymb = fftshift(ifft(X));
 
    % Transmitted signal is a sum of the delayed real, imag symbols
    symBuf = [symBuf(numFFT/2+1:end); complex(zeros(numFFT/2,1))];
    symBuf(KF+(1:KF)) = symBuf(KF+(1:KF)) + txSymb;
 
    % Compute power spectral density (PSD)
    currSym = complex(symBuf(1:KF));
    [specFBMC, fFBMC] = periodogram(currSym, hann(KF, 'periodic'), KF*2, 1);
    sumFBMCSpec = sumFBMCSpec + specFBMC;
 
    % Store transmitted signals for all symbols
    txSigAll(:,symIdx) = currSym;
end
 
%Scatterplot of modulated data (constellation diagram)
scatterplot(modData);title('modulated data')
 
% Plot power spectral density
sumFBMCSpec = sumFBMCSpec/mean(sumFBMCSpec(1+K+2*numGuards*K:end-2*numGuards*K-K));
plot(fFBMC-0.5,10*log10(sumFBMCSpec));
grid on
axis([-0.5 0.5 -180 10]);
xlabel('Normalized frequency');
ylabel('PSD (dBW/Hz)')
title(['FBMC, K = ' num2str(K) ' overlapped symbols'])
set(gcf, 'Position', figposition([15 50 30 30]));
 
% Process symbol-wise
BERR=zeros(size(snrdB));
for ii = 1:length(snrdB)
      BER = comm.ErrorRate;
for symIdx = 1:numSymbols
    rxSig = txSigAll(:, symIdx);
 
    % Add WGN
    rxNsig = awgn(rxSig, ii, 'measured');
 
    % Perform FFT
    rxf = fft(fftshift(rxNsig));
 
    % Matched filtering with prototype filter
    rxfmf = filter(Hk, 1, rxf);
    % Remove K-1 delay elements
    rxfmf = [rxfmf(K:end); zeros(K-1,1)];


    % Remove guards    
    rxfmfg = rxfmf(numGuards*K+1:end-numGuards*K);
 
    % OQAM post-processing
    %  Downsample by 2K, extract real and imaginary parts
    if rem(symIdx, 2)
        % Imaginary part is K samples after real one
        r1 = real(rxfmfg(1:2*K:end));
        r2 = imag(rxfmfg(K+1:2*K:end));
        rcomb = complex(r1, r2);
    else
        % Real part is K samples after imaginary one
        r1 = imag(rxfmfg(1:2*K:end));
        r2 = real(rxfmfg(K+1:2*K:end));
        rcomb = complex(r2, r1);
    end
    %  Normalize by the upsampling factor
    rcomb = (1/K)*rcomb;
 
    % De-mapper: Perform hard decision
    rxBits(:, symIdx) = qamdemod(rcomb, 2^bitsPerSubCarrier, ...
        'OutputType', 'bit', 'UnitAveragePower', true);
end
    % Measure BER with appropriate delay
    BER.ReceiveDelay = bitsPerSubCarrier*KL;
    ber = BER(inpData(:), rxBits(:));
    BERR(ii) = ber(1); %to keep BER for each SNR in variable BERR
end
 
%Scatterplot of received data (constellation diagram)
scatterplot(rcomb);title('received data')
 
%plot SER vs SNR in logarithm scale
Ser = BERR * bitsPerSubCarrier;
formula_Ser = 3/2*erfc(sqrt(0.1*(10.^(snrdB/10))));
figure
plot(snrdB,log10(Ser),'-x'); 		%numerics
hold on;
plot(snrdB,log10(formula_Ser),'o'); %analytics
legend('numerics','analytics');
xlabel('SNR, dB')
ylabel('Lg(SER)')
title('Symbol Error Rate (b)');
 
%plot BER vs SNR in logarithm scale
formula_Ber = (1/4.0)*3/2*erfc(sqrt(4*0.1*(10.^(snrdB/10))));
figure
plot(snrdB,log10(BERR),'-x');       %numerics
hold on;
plot(snrdB,log10(formula_Ber),'o'); %analytics
title('BER vs SNR (b)');
ylabel('Normalised BER');
xlabel('SNR (dB)');
grid on
legend('numerics','analytics');
 
%eyediagrams of transmitter and receiver
eyediagram(modData,num_samples_each_trace) %transmitter
eyediagram(rcomb,num_samples_each_trace) 	 %receiver
 	
%Q-factor
Q = 20*log10(sqrt(2)*(erfcinv(2*BERR)));
formula_Q = 20*log10(sqrt(2)*(erfcinv(2*formula_Ber)));
figure
plot(snrdB,Q,'-x');        %numerics
hold on;
plot(snrdB,formula_Q,'o'); %analytics
title('Q factor vs SNR (b)');
ylabel('Q factor');
xlabel('SNR (dB)');
grid on
legend('numerics','analytics');
 
%BW efficiency
S=ceil(length(modData)/(numSymbols/bitsPerSubCarrier));    
BW = S/(S+K-(1/2));  %BW efficiency formula of FBMC
fprintf('Bandwidth efficiency is %.2f bit/s/Hz\n',BW)
 
% Restore RNG state
rng(s);

%% OFDM Power Spectral Density
for symIdx = 1:numSymbols

    inpData2 = randi([0 1], bitsPerSubCarrier*L, 1);
    modData = qammod(inpData2, 2^bitsPerSubCarrier, ...
        'InputType', 'Bit', 'UnitAveragePower', true);

    symOFDM = [zeros(numGuards,1); modData; zeros(numGuards,1)];
    ifftOut = sqrt(numFFT).*ifft(ifftshift(symOFDM));

    [specOFDM,fOFDM] = periodogram(ifftOut, rectwin(length(ifftOut)), ...
        numFFT*2, 1, 'centered');
    sumOFDMSpec = sumOFDMSpec + specOFDM;
end

% Plot power spectral density (PSD) over all subcarriers
sumOFDMSpec = sumOFDMSpec/mean(sumOFDMSpec(1+2*numGuards:end-2*numGuards));
figure;
plot(fOFDM,10*log10(sumOFDMSpec));
grid on
axis([-0.5 0.5 -180 10]);
xlabel('Normalized frequency');
ylabel('PSD (dBW/Hz)')
title(['OFDM, numFFT = ' num2str(numFFT)])
set(gcf, 'Position', figposition([46 50 30 30]));

toc
