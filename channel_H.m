% This function generates the channel Matrix Coefficients
% channel_H(i, j, t) = hij during the time interval 't' in which
%                      a frame is transmitted
%                      Thus, this is for a slow fading channel
%                      (frame static channel coefficients)

function [v] = channel_H(i, j, t, itr, N, fD);


phimax=2*pi*2/360;
alphaBS=2*pi*90/360;
alphaMS=2*pi*90/360;
alphav=2*pi*180/360;
SNR=10^(17/10);


if (itr == 1)
    deltaBS=1/2;
    deltaMS=1/2;
end

if (itr == 2)
        deltaBS=1;
        deltaMS=1;
end

if (itr == 3)
    deltaBS=30;
    deltaMS=3;
end

teta=2*pi*rand(1,N);


% kappa = 30; mu = 180;
phi=[
   3.03470037473462
  -3.20686453918495
   2.74685103894036
  -2.97465139867191
  -3.05510972018990
  -3.21044033537730
  -2.89550049059082
   2.91232099977409
  -2.74973414161808
  -3.05399746121793
   2.90144584407880
  -3.03999800647819
  -3.14555442569988
  -2.91462440384403
  -3.13965293772270
  -3.13102673664680
  -3.01283924085586
   3.03155817522709
  -3.32455978845568
  -3.24487162849587
  -3.20002334529738
   3.16982714092315
  -3.06017006446165
  -2.89806096831352
   2.90093813889891
];



%t=32;
sum=0;

if i==1
    if j==1
        for n=1:N
            
            an=exp(sqrt(-1)*pi*deltaBS*(cos(alphaBS)+phimax*sin(alphaBS)*sin(phi(n))));
            bn=exp(sqrt(-1)*pi*deltaMS*cos(phi(n)-alphaMS));
            fn=fD*cos(phi(n)-alphav);
            sum=sum+an*bn*exp(sqrt(-1)*(2*pi*fn*t+teta(n)));
        end
        v=sum/sqrt(N);
    elseif j==2
        for n=1:N
            
            an=exp(sqrt(-1)*pi*deltaBS*(cos(alphaBS)+phimax*sin(alphaBS)*sin(phi(n))));
            bn=exp(sqrt(-1)*pi*deltaMS*cos(phi(n)-alphaMS));
            fn=fD*cos(phi(n)-alphav);
            sum=sum+conj(an)*bn*exp(sqrt(-1)*(2*pi*fn*t+teta(n)));
        end
        v=sum/sqrt(N);
    end    
elseif i==2
    if j==1
        for n=1:N
            
            an=exp(sqrt(-1)*pi*deltaBS*(cos(alphaBS)+phimax*sin(alphaBS)*sin(phi(n))));
            bn=exp(sqrt(-1)*pi*deltaMS*cos(phi(n)-alphaMS));
            fn=fD*cos(phi(n)-alphav);
            sum=sum+an*conj(bn)*exp(sqrt(-1)*(2*pi*fn*t+teta(n)));
        end 
        v=sum/sqrt(N);
    elseif j==2
        for n=1:N
            
            an=exp(sqrt(-1)*pi*deltaBS*(cos(alphaBS)+phimax*sin(alphaBS)*sin(phi(n))));
            bn=exp(sqrt(-1)*pi*deltaMS*cos(phi(n)-alphaMS));
            fn=fD*cos(phi(n)-alphav);
            sum=sum+conj(an)*conj(bn)*exp(sqrt(-1)*(2*pi*fn*t+teta(n)));
        end     
        v=sum/sqrt(N);
    end
end   
