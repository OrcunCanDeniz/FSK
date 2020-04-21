%Orcun Can Deniz 
%B/M-FSK Simulation Script 

%% Modulate, Send, Demodulate B/M FSK Signal
clc;clear all;close all
freqs_for_visualization = true;
M=8;%change this value to have FSKs with different order
sym_num = 10; %number of symbols in stream
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

colors = rand(M,3);
fc = 20e+3;
freqs=[];

if freqs_for_visualization == false 
    symDur=0.0005;
    freq_increment = 1/(2*symDur);
    freqs= linspace(fc,fc+M*freq_increment,M)
else
    symDur = 0.75
    freqs = linspace(M,M*M,M);
end
fs=freqs(end)*4; %set Fsampling 
Ts=1/fs;
mod = [];
figure;
vis_data = [data data(end)];%preprocess data stream for visualization
subplot(2,1,1);stairs(vis_data, 'LineWidth', 3);xlim([1 numel(vis_data)]);title("Data Stream")

for i=1:length(data)
    tSym = (i-1)*symDur:Ts:(symDur*i);% symbol duration
    tmp = cos(2*pi*tSym* freqs(data(i)+1));%%modulate data with corresponding frequency
    mod = [mod tmp];%append the symbol duration to complete signal
    subplot(2,1,2);plot(tSym,tmp,'color', colors(data(i)+1,:));
    title("Modulated Signal @TX");xlabel('Time');ylabel("Amplitude[V]")
    hold on; axis([0 tSym(end) -2 2]);
    hold on
end

F = linspace(-fs/2,fs/2, numel(mod));%frequency range for plotting
t = linspace(0,numel(data),numel(mod));%time range for plotting
fft_TX = fftshift(abs(fft(mod)));
fft_TX = 1/numel(fft_TX)*fft_TX;%Mag spectrum of signal
figure;
plot(F,fft_TX);title("Mag. Spectrum of Signal @TX");xlabel('Freq. [Hz]');ylabel("Magnitude[V]")


figure;
noisy_sig = awgn(mod,-2,'measured');%pass signal through awgn ch.
subplot(2,1,1);plot(t,noisy_sig);title("FSK Signal @RX");xlabel('Time');ylabel("Amplitude[V]")

fft_NO = fftshift(abs(fft(noisy_sig)));
fft_NO = 1/numel(fft_NO)*fft_NO;%spectrum @ RX
subplot(2,1,2);plot(F,fft_NO);title("Mag. Spectrum of Signal @RX");
xlabel('Freq. [Hz]');ylabel("Magnitude[V]")

n_sym=numel(data);
n_mod=numel(mod); 
sPerSymb = n_mod/n_sym;%sample per symbol
n_freqs=numel(freqs)
correlator = zeros(n_sym,n_freqs);%preallocation for demoding
%each row is repsenting period for one symbol
%columns are representing the correlator output for a freq in a sym period

for i=1:n_sym
    tSym = (i-1)*symDur:Ts:(symDur*i);
    for j=1:n_freqs
        out = noisy_sig((i-1)*sPerSymb+1:i*sPerSymb).*cos(2*pi*tSym* freqs(j));
        %multiplying with cos waves at expected freqs in above line
        correlator(i,j) = trapz(tSym,out);
        %integrating multiplied signal and inserting to a matrix
    end
end

[Maxes,I] = max(correlator,[],2);
%in each row, max valued element's index-1 corresponds to symbol
decoded_data = (I-1)'
data


figure;
vis_dec_data = [decoded_data decoded_data(end)];
stairs(vis_dec_data, 'LineWidth', 2,'color',[0 1 0]);xlim([1 numel(vis_dec_data)]);title("Decoded Data Stream")

%% Compare BER and Power Efficiencies
EbNo = (-5:10)';

ber4 = berawgn(EbNo,'fsk',4,'coherent');
ber8 = berawgn(EbNo,'fsk',8,'coherent');
figure;
semilogy(EbNo,[ber4 ber8], 'LineWidth',2); xlabel("dB"); ylabel("BER");
legend('4FSK','8FSK') ;

%% Compare BW Efficiencies
Ms = (2:256);
bwEff = 2.*log2(Ms)./Ms
figure;
semilogy(Ms,[bwEff], 'LineWidth',4); xlabel("M"),
title("MFSK Bandwith Efficiency");ylabel("Bandwidth Effieciency [b/s/Hz]");xlim([2 256])