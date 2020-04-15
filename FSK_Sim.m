%Orcun Can Deniz 16014103
%B/M-FSK Simulation Script as Comm2 Assignment
clc;clear all;close all

M=8;%change this value to have different order FSKs
switch M
    case 2
        data=[0,1,0,1,0,1,0,1];
    case 4
        data=[0,3,1,2,0,3,1,2];
    case 8
        data=[7,0,6,1,5,2,4,3];
    otherwise
        data=randi([0 M-1],8,1)'
end
load("colors.mat")
freqs = linspace(M,M*6,M); %Set frequencies
fs=freqs(end)*4; %set Fsampling 
Ts=1/fs;
mod = [];
figure;
vis_data = [data data(end)];%preprocess data stream for visualization
subplot(2,1,1);stairs(vis_data, 'LineWidth', 2);xlim([1 numel(vis_data)]);title("Data Stream")

for i=1:length(data)
    tSym = i-1:Ts:i-Ts; % symbol duration
    tmp = cos(2*pi*tSym* freqs(data(i)+1));%get the wave at the frequency of the symbol 
    mod = [mod tmp];%append the symbol duration to complete signal
    subplot(2,1,2);plot(tSym,tmp,'color', colors(data(i)+1,:));
    title("Modulated Signal @TX");xlabel('Time');ylabel("Amplitude[V]")
    hold on; axis([0 numel(data) -2 2]);
    hold on;
end

F = linspace(-fs/10,fs/10, numel(mod));%frequency range for plotting
t = linspace(0,numel(data),numel(mod));%time range for plotting
fft_TX = fftshift(abs(fft(mod)));
fft_TX = 1/numel(fft_TX)*fft_TX;%Mag spectrum of signal
figure;
plot(F,fft_TX);title("Mag. Spectrum of Signal @RX");xlabel('Freq. [Hz]');ylabel("Magnitude[V]")


figure;
noisy_sig = awgn(mod,-5,'measured');%pass signal through awgn ch.
subplot(2,1,1);plot(t,noisy_sig);title("FSK Signal @RX");xlabel('Time');ylabel("Amplitude[V]")

fft_NO = fftshift(abs(fft(noisy_sig)));
fft_NO = 1/numel(fft_NO)*fft_NO;%spectrum @ RX
subplot(2,1,2);plot(F,fft_NO);title("Mag. Spectrum of Signal @RX");
xlabel('Freq. [Hz]');ylabel("Magnitude[V]")

n_sym=numel(data);
n_freqs=numel(freqs)
correlator = zeros(n_sym,n_freqs);%preallocation for demoding
%each row is repsenting period for one symbol
%columns are representing the correlator output for a freq in a sym period

for i=1:n_sym
    tSym = i-1:Ts:i-Ts;
    for j=1:n_freqs
        out = noisy_sig((i-1)*fs+1:i*fs).*cos(2*pi*tSym* freqs(j));
        %multiplying with cos waves at expected freqs in above line
        correlator(i,j) = trapz(tSym,out);
        %integrating multiplied signal and inserting to a matrix
    end
end

[M,I] = max(correlator,[],2);
%in each row, max valued element's index-1 corresponds to symbol
decoded_data = (I-1)'
data

figure;
vis_dec_data = [decoded_data decoded_data(end)];
stairs(vis_dec_data, 'LineWidth', 2,'color',[0 1 0]);xlim([1 numel(vis_dec_data)]);title("Decoded Data Stream")

EbNo = (-5:10)';

ber4 = berawgn(EbNo,'fsk',4,'coherent');
ber8 = berawgn(EbNo,'fsk',8,'coherent');
figure;
semilogy(EbNo,[ber4 ber8], 'LineWidth',2); xlabel("dB"); ylabel("BER");
legend('4FSK','8FSK') ;
