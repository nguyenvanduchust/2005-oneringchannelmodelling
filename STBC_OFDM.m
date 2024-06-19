clear all

% STBC OFDM in frequency selective channel based on one ring model
diary text_itr4_fd_60.txt

load channel_profile.am;

rho = channel_profile;

R = 150; % tau_max = (2*R/c)=1000 ns
%**********************************
% itr

% itr= 1; 1. run: delta_BS/lambda= 1/2  and delta_MS/lambda= 1/2  
% itr= 2; 2. run: delta_BS/lambda= 1  and delta_MS/lambda= 1
% itr= 3; 3. run: delta_BS/lambda= 30 and delta_MS/lambda=3
% itr= 4; 4. run: delta_BS/lambda= 30 and delta_MS/lambda= 1/2
% itr= 5; 5. run: delta_BS/lambda= 1/2 and delta_MS/lambda=3
% itr= 6; 6. run: delta_BS/lambda= 3 and delta_MS/lambda= 1

%**********************************


NFFT = 64;                  %FFT length
G = 21;                  

M_ary =16;                  %Multilevel of M-ary symbol

 

fD = 60.0;

N_P  = length(rho);             %Length of the CIR


t_a = 50*10^(-9);               %Sampling duration of HiperLAN/2

symbol_duration = NFFT * t_a;   %OFDM symbol duration

itr = 4;
N = 80;

mse = [];
NofOFDMSymbol = 30;            %Number of data and pilot OFDM symbol

teta=2*pi*rand(1,N);

length_data =  NofOFDMSymbol * NFFT;  % the total data length

y_16QAM=[1+i;-1+i;1-i;-1-i; 3+i;-3+i;3-i; -3-i; 1+3*i; -1+3*i; 1-3i; -1-3*i; 3+3*i; -3+3*i; 3-3*i; -3-3*i ];
%-------------
% Source bites
%-------------

source_data1 = randint(length_data,4);
source_data2 = randint(length_data,4);
%--------------------
% bit to symbol coder
%--------------------

symbols1 = bi2de(source_data1);  
symbols2 = bi2de(source_data2);
%----------------------------
% 16 QAM modulator in base band
%----------------------------

QAM_Symbol1=qam16(symbols1);
QAM_Symbol2=qam16(symbols2);
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
    
    clear QAM_tem;

end;


%----------------------------
% STC coding
%--------------------------------
stc_n1=[];
stc_cn1=[];
sct_n2=[];
stc_cn2=[];
stc_n1=Data_Pattern1;
stc_cn1=-conj(Data_Pattern2);
stc_n2=Data_Pattern2;
stc_cn2=conj(Data_Pattern1);


clear Data_Pattern1 Data_Pattern2;


%--------------------------------------------------------------------------
% Transmitted signal of trasmitt antena 1 to receive antenna 1 (channel: h11)
% and transmitted signal of trasmitt antena 1 to receive antenna 2 (channel: h12)
%--------------------------------------------------------------------------

SER_R = [];


snr_min =17;
snr_max =19.0;
step = 2.0;
SER=[];
% snr=10;
 for snr = snr_min:step:snr_max; 
       
     SER_t=0;
    
     snr = snr - 10*log10((NFFT-G)/NFFT); %Miss matching effect
  
     if snr<11.73 % to make sure the length of symbol is long enough and SER is correct
         NT=100
     else
         NT=100000
     end
     
    for  t_i=1:NT %(1000*1920)

rs11_frame1 = [];
rs11_frame2 = [];
rs12_frame1 = [];
rs12_frame2 = [];

h11_frame = [];
h12_frame = [];

t = 0.1;

for i=0:NofOFDMSymbol-1;
    
   OFDM_signal_tem11 = OFDM_Modulator(stc_n1(i+1,:),NFFT,G); % assuming that two adjacent symbols experience the same fading
   
   OFDM_signal_tem12 = OFDM_Modulator(stc_cn1(i+1,:),NFFT,G);
   
   % OFDM signal from the first antenna is created
   [h11] = OneRingSelectiveCh(1, 1, t, itr, N, fD, t_a, R, teta, rho);
   
   h11_frame = [h11_frame; h11];
   
   
   rs11_1 = conv(OFDM_signal_tem11, h11);
   rs11_2 = conv(OFDM_signal_tem12, h11); % the second symbol
   % The received signal over multhipath channel is created
   rs11_1 = awgn(rs11_1,snr,'measured','dB');
   rs11_2 = awgn(rs11_2,snr,'measured','dB');
   
   rs11_frame1 = [rs11_frame1; rs11_1];
  
   rs11_frame2 = [rs11_frame2; rs11_2];
   
%-----------------------------------------------------------------------
% Transmitted signal of antenna 1 to receive antenna 2 (channel: h12)
%-----------------------------------------------------------------------



   [h12] = OneRingSelectiveCh(1, 2, t, itr, N, fD, t_a, R, teta, rho);
   h12_frame = [h12_frame; h12];
   
   
   rs12_1 = conv(OFDM_signal_tem11, h12);
   rs12_2 = conv(OFDM_signal_tem12, h12);
   
   % The received signal over multhipath channel is created
   
   rs12_1 = awgn(rs12_1,snr,'measured','dB');
   rs12_2 = awgn(rs12_2,snr,'measured','dB');
   
   % The received signal over multhipath channel with additive noise is created
   
   rs12_frame1 = [rs12_frame1; rs12_1];
   rs12_frame2 = [rs12_frame2; rs12_2];
   
   t = t + symbol_duration;
    clear OFDM_signal_tem11, OFDM_signal_tem12;
     
end;

%------------------------------------------------------------------------
% Transmitted signal of antenna 2 to receive antenna 1 (channel: h21)
% and transmitted signal of antenna 2 to receive antenna 2 (channel: h22)
%-----------------------------------------------------------------------

t = 0.1;
rs21_frame1 = [];
rs21_frame2 = [];
rs22_frame1 = [];
rs22_frame2 = [];

h21_frame = [];
h22_frame = [];

for i=0:NofOFDMSymbol-1;
    
   OFDM_signal_tem21 = OFDM_Modulator(stc_n2(i+1,:),NFFT,G);
   
   OFDM_signal_tem22 = OFDM_Modulator(stc_cn2(i+1,:),NFFT,G);
   
   % OFDM signal from the second antenna is created
   [h21] = OneRingSelectiveCh(2, 1, t, itr, N, fD, t_a, R, teta, rho);
   h21_frame = [h21_frame; h21];
   
   
   rs21_1 = conv(OFDM_signal_tem21, h21);
   rs21_2 = conv(OFDM_signal_tem22, h21);
   % The received signal over multhipath channel is created
   rs21_1 = awgn(rs21_1,snr,'measured','dB');
   
   rs21_2 = awgn(rs21_2,snr,'measured','dB');
   
   rs21_frame1 = [rs21_frame1; rs21_1];
   rs21_frame2 = [rs21_frame2; rs21_2];
   
      
%------------------------------------------------------------------------
% Transmitted signal of antenna 2 to receive antenna 2 (channel: h22)
%-----------------------------------------------------------------------


  
   % OFDM signal from the second antenna is created

   [h22] = OneRingSelectiveCh(2, 2, t, itr, N, fD, t_a, R, teta, rho);
   
   h22_frame = [h22_frame; h22];
   
   
   rs22_1 = conv(OFDM_signal_tem21, h22);
   rs22_2 = conv(OFDM_signal_tem22, h22);  
   
   % The received signal over multhipath channel is created
   
   rs22_1 = awgn(rs22_1,snr,'measured','dB');
   rs22_2 = awgn(rs22_2,snr,'measured','dB');
   
   % The received signal over multhipath channel with additive noise is created
   rs22_frame1 = [rs22_frame1; rs22_1];
   rs22_frame2 = [rs22_frame2; rs22_2];
   
   t = t + symbol_duration;
    
   clear OFDM_signal_tem21,OFDM_signal_tem22;
   
end;

%-----------------
% Recever 1: 
%-----------------

Receiver_Data_n1 = []; %  Prepare a matrix for reveived data symbols
Receiver_Data_n2 = []; % the second symbol

SignalPostFFT11 = [];
SignalPostFFT12 = [];
SignalPostFFT21 = [];
SignalPostFFT22 = [];

d11 = []; % Received signal of the first receive antenna
d12 = [];
d21 = []; % Received signal of the second receive antenna
d22 = []; 
data_symbol_1 = [];
data_symbol_2 = [];
s1_est=[];
s2_est=[];

for i=1:NofOFDMSymbol;
   
          rs1_1i = rs11_frame1(i,:) + rs21_frame1(i,:); % first symbol
          rs1_2i = rs12_frame1(i,:) + rs22_frame1(i,:);

          rs2_1i = rs11_frame2(i,:) + rs21_frame2(i,:);  %second symbol
          rs2_2i = rs12_frame2(i,:) + rs22_frame2(i,:);
          
        %-----------------------------------------------------------
        % OFDM demodulator
        %-----------------------------------------------------------
        
        % the first received symbol
        
        Demodulated_signal1_1i = OFDM_Demodulator(rs1_1i,NFFT,NFFT,G);  
        SignalPostFFT11 = [SignalPostFFT11; Demodulated_signal1_1i];
        
        Demodulated_signal1_2i = OFDM_Demodulator(rs1_2i,NFFT,NFFT,G);  
        SignalPostFFT12 = [SignalPostFFT12; Demodulated_signal1_2i];
        %-------------------------------------------------
        % the second received symbol
        
        Demodulated_signal2_1i = OFDM_Demodulator(rs2_1i,NFFT,NFFT,G); 
        SignalPostFFT21 = [SignalPostFFT21; Demodulated_signal2_1i];
        
        Demodulated_signal2_2i = OFDM_Demodulator(rs2_2i,NFFT,NFFT,G); 
        SignalPostFFT22 = [SignalPostFFT22; Demodulated_signal2_2i];
       
        
        h11_i = h11_frame(i,:);
        H11_i = fft([h11_i,zeros(1,NFFT-N_P)]);
        h12_i = h12_frame(i,:);
        H12_i = fft([h12_i,zeros(1,NFFT-N_P)]);
        
        h21_i = h21_frame(i,:);
        H21_i = fft([h21_i,zeros(1,NFFT-N_P)]);
        h22_i = h22_frame(i,:);
        H22_i = fft([h22_i,zeros(1,NFFT-N_P)]);
        
  %----------------------------
  % STC decoding
  %---------------------------   
  r0=[];
  r1=[];
  r2=[];
  r3=[];
  
    r0=Demodulated_signal1_1i;
    r1=Demodulated_signal2_1i;
    r2=Demodulated_signal1_2i;
    r3=Demodulated_signal2_2i;
    
    s1_esti= conj(H11_i).*r0+H21_i.*conj(r1)+conj(H12_i).*r2+H22_i.*conj(r3);
    s2_esti= conj(H21_i).*r0-H11_i.*conj(r1)+conj(H22_i).*r2-H12_i.*conj(r3);
  
       clear r0;clear r1; clear r2; clear r3;    
 

%    
 HK=(conj(H11_i).*H11_i+ H21_i.*conj(H21_i)+conj(H12_i).*H12_i+H22_i.*conj(H22_i));

 for  i=1:length(s1_esti)
     
        for ll=1:length(y_16QAM)
              d(ll)=(HK(i)-1)*(abs(y_16QAM(ll))).^2+(abs(s1_esti(i)-y_16QAM(ll))).^2;
        end
          [D_min index]=min(d);
          s1_esti(i)=y_16QAM(index);
 end
      clear index
for  i=1:length(s2_esti)
    
        for ll=1:length(y_16QAM)
              d(ll)=(abs(HK(i))-1)*(abs(y_16QAM(ll))).^2+(abs(s2_esti(i)-y_16QAM(ll))).^2;
        end
          [D_min index]=min(d);
          s2_esti(i)=y_16QAM(index);
 end
 
    s1_est=[s1_est s1_esti];
    s2_est=[s2_est s2_esti];
    
        
   end  % end of line 266

        s1 = ddemodce(s1_est,1,1,'qam',M_ary);
        s2 = ddemodce(s2_est,1,1,'qam',M_ary);
       
        data_symbol_1 = s1';
        data_symbol_2 = s2';
        [number1, ratio1] = symerr(symbols1,data_symbol_1);
        [number2, ratio2] = symerr(symbols2,data_symbol_2);
        
        SER_R = (ratio1+ratio2)/2; 
        
        SER_t=SER_t+SER_R;
        
 
    end
       tE=SER_t/t_i;
       SER = [SER, tE]; 
      %SER=tE;
      
  end  % end of snr

  

 psnr = snr_min-3:step:snr_max-3; %for two antennas divide the total energy
        
    
semilogy(psnr, SER,'ro-');


xlabel('SNR in dB');

ylabel('SER');



hold on