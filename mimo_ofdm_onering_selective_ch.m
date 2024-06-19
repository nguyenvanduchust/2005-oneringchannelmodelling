%--------------------------------------------------------------------------

%-------------------------------------------------------------------------

%clc;
clear all;

load channel_profile.am;

rho = channel_profile;



R = 150; % tau_max = (2*R/c)=1000 ns

NFFT = 64;                  %FFT length
G = 21;                  

M_ary =16;                  %Multilevel of M-ary symbol

 

fD = 91.0;





N_P  = length(rho);             %Length of the CIR


t_a = 50*10^(-9);               %Sampling duration of HiperLAN/2

symbol_duration = NFFT * t_a;   %OFDM symbol duration

itr =3;
N = 80;

NofOFDMSymbol = 30;            %Number of data and pilot OFDM symbol
mse = [];


teta=2*pi*rand(1,N);


length_data =  NofOFDMSymbol * NFFT;  % the total data length

%-------------
% Source bites
%-------------
source_data1 = randi([0 1],length_data,4);
source_data2 = randi([0 1],length_data,4);
%--------------------
% bit to symbol coder
%--------------------

symbols1 = bi2de(source_data1);  
symbols2 = bi2de(source_data2);  
%----------------------------
% 16 QAM modulator in base band
%----------------------------------------------
%QAM_Symbol1 = dmodce(symbols1,1,1,'qam',M_ary);
%QAM_Symbol2 = dmodce(symbols2,1,1,'qam',M_ary);
QAM_Symbol1=qammod(symbols1,M_ary);
QAM_Symbol2=qammod(symbols2,M_ary);

%-----------------------------------------------
% Preparing data pattern
%
%-----------------------------------------------

Data_Pattern1 = []; % Transmitted Signal before IFFT
m = 0;
for i=0:NofOFDMSymbol-1;
    QAM_tem = [];
    for n=1:NFFT;
          QAM_tem = [QAM_tem,QAM_Symbol1(i*NFFT+n)];
    end;
    Data_Pattern1 = [Data_Pattern1;QAM_tem];
    
    clear QAM_tem;

end;

Data_Pattern2 = []; % Transmitted Signal before IFFT
m = 0;
for i=0:NofOFDMSymbol-1;
    QAM_tem = [];
    for n=1:NFFT;
          QAM_tem = [QAM_tem,QAM_Symbol2(i*NFFT+n)];
    end;
    Data_Pattern2 = [Data_Pattern2;QAM_tem];
    
    %clear QAM_tem;

end;




%--------------------------------------------------------------------------
% Transmitted signal of trasmitt antena 1 to receive antenna 1 (channel: h11)
% and transmitted signal of trasmitt antena 1 to receive antenna 2 (channel: h12)
%--------------------------------------------------------------------------

SER_R = [];


snr_min =3.0;
snr_max =39.0;
step = 4.0;
for snr = snr_min:step:snr_max; 
     SER_t=0;

     for  t_i=1:100 %(1000*1920)
snr = snr - 10*log10((NFFT-G)/NFFT); %Miss matching effect

rs11_frame = [];
rs12_frame = [];


h11_frame = [];
h12_frame = [];

t = 0.1;

for i=0:NofOFDMSymbol-1;
   OFDM_signal_tem1 = OFDM_Modulator(Data_Pattern1(i+1,:),NFFT,G);
   
   % OFDM signal from the first antenna is created
   [h11] = OneRingSelectiveCh(1, 1, t, itr, N, fD, t_a, R, teta, rho);
   
   h11_frame = [h11_frame; h11];
   
   
   rs11 = conv(OFDM_signal_tem1, h11);
   % The received signal over multhipath channel is created
   rs11 = awgn(rs11,snr,'measured','dB');
   
   
   rs11_frame = [rs11_frame; rs11];
  
   
   
%-----------------------------------------------------------------------
% Transmitted signal of antenna 1 to receive antenna 2 (channel: h12)
%-----------------------------------------------------------------------



   [h12] = OneRingSelectiveCh(1, 2, t, itr, N, fD, t_a, R, teta, rho);
   h12_frame = [h12_frame; h12];
   
   
   rs12 = conv(OFDM_signal_tem1, h12);
   % The received signal over multhipath channel is created
   rs12 = awgn(rs12,snr,'measured','dB');
   
   % The received signal over multhipath channel with additive noise is created
   rs12_frame = [rs12_frame; rs12];
   
   t = t + symbol_duration;
    clear OFDM_signal_tem1;
   
   
end;


%------------------------------------------------------------------------
% Transmitted signal of antenna 2 to receive antenna 1 (channel: h21)
% and transmitted signal of antenna 2 to receive antenna 2 (channel: h22)
%-----------------------------------------------------------------------

t = 0.1;
rs21_frame = [];
rs22_frame = [];

h21_frame = [];
h22_frame = [];

for i=0:NofOFDMSymbol-1;
   OFDM_signal_tem2 = OFDM_Modulator(Data_Pattern2(i+1,:),NFFT,G);
   
   
   % OFDM signal from the second antenna is created
   [h21] = OneRingSelectiveCh(2, 1, t, itr, N, fD, t_a, R, teta, rho);
   h21_frame = [h21_frame; h21];
   
   
   rs21 = conv(OFDM_signal_tem2, h21);
   % The received signal over multhipath channel is created
   rs21 = awgn(rs21,snr,'measured','dB');
   
   
   rs21_frame = [rs21_frame; rs21];
   
   
   
%------------------------------------------------------------------------
% Transmitted signal of antenna 2 to receive antenna 2 (channel: h22)
%-----------------------------------------------------------------------


  
   % OFDM signal from the second antenna is created

   [h22] = OneRingSelectiveCh(2, 2, t, itr, N, fD, t_a, R, teta, rho);
   
   h22_frame = [h22_frame; h22];
   
   
   rs22 = conv(OFDM_signal_tem2, h22);
   % The received signal over multhipath channel is created
   rs22 = awgn(rs22,snr,'measured','dB');
   
   % The received signal over multhipath channel with additive noise is created
   rs22_frame = [rs22_frame; rs22];
   
   t = t + symbol_duration;
    
   clear OFDM_signal_tem2;
   
end;





%-----------------
% Recever 1: 
%-----------------

Receiver_Data = []; %  Prepare a matrix for reveived data symbols

SignalPostFFT1 = [];
SignalPostFFT2 = [];
SER_t=0;

d1 = []; % Received signal of the first receive antenna
d2 = []; % Received signal of the second receive antenna
data_symbol_1 = [];
data_symbol_2 = [];

for i=1:NofOFDMSymbol;
   
          rs1_i = rs11_frame(i,:) + rs21_frame(i,:);
          rs2_i = rs12_frame(i,:) + rs22_frame(i,:);

        
        %-----------------------------------------------------------
        % OFDM demodulator
        %-----------------------------------------------------------
        Demodulated_signal1_i = OFDM_Demodulator(rs1_i,NFFT,NFFT,G);  
        SignalPostFFT1 = [SignalPostFFT1; Demodulated_signal1_i];

       %--------------------------------------------------------------
        Demodulated_signal2_i = OFDM_Demodulator(rs2_i,NFFT,NFFT,G); 
        SignalPostFFT2 = [SignalPostFFT2; Demodulated_signal2_i];
        
       
        
        h11_i = h11_frame(i,:);
        H11_i = fft([h11_i,zeros(1,NFFT-N_P)]);
        h12_i = h12_frame(i,:);
        H12_i = fft([h12_i,zeros(1,NFFT-N_P)]);
        
        h21_i = h21_frame(i,:);
        H21_i = fft([h21_i,zeros(1,NFFT-N_P)]);
        h22_i = h22_frame(i,:);
        H22_i = fft([h22_i,zeros(1,NFFT-N_P)]);
   
        
        
        
        
       
        d1_i = [];
        d2_i = [];
        for k = 1:NFFT;
            H_k = [H11_i(k),H21_i(k); H12_i(k),H22_i(k)];
            y = [Demodulated_signal1_i(k);Demodulated_signal2_i(k)];
            x = inv(H_k) * y;
            d1_i = [d1_i,x(1)];
            d2_i = [d2_i,x(2)];
        end;
        
        %--------------------------
        % Demodulated signal 
        %--------------------------
        d1 = [d1; d1_i];
        
        d2 = [d2; d2_i];
        
        %demodulated_symbol_1i = ddemodce(d1_i,1,1,'qam',M_ary);
        %demodulated_symbol_2i = ddemodce(d2_i,1,1,'qam',M_ary);
        demodulated_symbol_1i = qamdemod(d1_i,M_ary);
        demodulated_symbol_2i = qamdemod(d2_i,M_ary);

        data_symbol_1 = [data_symbol_1, demodulated_symbol_1i];
        
        data_symbol_2 = [data_symbol_2, demodulated_symbol_2i];
        

        
    end;
    
        data_symbol_1 = data_symbol_1';
        data_symbol_2 = data_symbol_2';
        [number1, ratio1] = symerr(symbols1,data_symbol_1);
        [number2, ratio2] = symerr(symbols2,data_symbol_2);
        
        %SER_R = [SER_R, (ratio1+ratio2)/2]; 
         tem = (ratio1+ratio2)/2; 
         SER_t=SER_t+tem;
       % SER_t=SER_t+SER_R;
        
    end
    
        
       tE=SER_t/t_i;
       SER_R = [SER_R, tE]; 
       %SER=tE;

end;

    
snr = snr_min:step:snr_max;
        
    

 figure(1);
 
 plot(real(d1),imag(d1),'b+');

 
 axis([-5 5 -5 5]);
 
 
xlabel('In-phase component of interference cancelled symbols');

 ylabel('Quadrature component of interference cancelled symbols');
 
 
 figure(2);

semilogy(snr, SER_R,'ro-');



xlabel('SNR in dB');

ylabel('SER');




