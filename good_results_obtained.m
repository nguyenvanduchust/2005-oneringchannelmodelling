clear all;



load ser_vd_NoOFDMSymbol20_NT_100_itr_4_firsttest.mat;
semilogy(psnr, SER,'ro-');
hold on;


load ser_vd_NoOFDMSymbol20_NT_100_itr_1_firsttest.mat;
semilogy(psnr, SER,'b*-');
%load ser_vd_NoOFDMSymbol20_NT_100_itr_1_secondtest.mat;
%semilogy(psnr, SER,'bs:');


load ser_vd_NoOFDMSymbol20_NT_100_itr_5_firsttest.mat;
semilogy(psnr, SER,'kv-');
load ser_vd_NoOFDMSymbol20_NT_100_itr_3_firsttest.mat;
semilogy(psnr, SER,'gs-');

xlabel('SNR in dB');

ylabel('SER');

legend('deltaBS=30, deltaMS=1/2 (itr = 4)','deltaBS=1/2, deltaMS=1/2 (itr = 1)',...
     'deltaBS=1/2, deltaMS=3 (itr = 5)','deltaBS=30 , deltaMS=3 (itr = 3)');
title('f_D = 60.0 Hz');
hold off;
grid on;

