 function [QAM16]=qam16(x)
% clear
% x=[1 5 6 9 8 4 11 12]
y=[1+i;-1+i;1-i;-1-i; 3+i;-3+i;3-i; -3-i; 1+3*i; -1+3*i; 1-3i; -1-3*i; 3+3*i; -3+3*i; 3-3*i; -3-3*i ];

for l=1:length(x);
      t=x(l)+1;

      QAM16(l,1)=y(t,1);

 end