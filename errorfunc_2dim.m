function [E]=errorfunc_2dim(phi)
header;
deltaBS=linspace(0,30,6);
deltaMS=linspace(0,3,6);
phi2=linspace(-pi,pi,100);
for i=1:length(deltaBS)
    for j=1:length(deltaMS)
        a=exp(sqrt(-1).*pi.*deltaBS(i).*(cos(alphaBS)+phimax.*sin(alphaBS).*sin(phi2)));
        b=exp(sqrt(-1).*pi.*deltaMS(j).*cos(phi2-alphaMS));
        rr=trapz(phi2,a.^2.*b.^2.*Mises(phi2,kappa,mu));
        sum=0;
        for n=1:N
            aa=exp(sqrt(-1).*pi.*deltaBS(i).*(cos(alphaBS)+phimax.*sin(alphaBS).*sin(phi(n))));
            bb=exp(sqrt(-1).*pi.*deltaMS(j).*cos(phi(n)-alphaMS));
            sum=sum+aa.^2.*bb.^2;
        end
        ss=sum/N;
        int1(j)=(abs(rr-ss))^p;    
    end
    int2(i)=trapz(deltaMS,int1);
end
deltaBS_max=deltaBS(length(deltaBS));
deltaMS_max=deltaMS(length(deltaMS));
E=(1/(deltaBS_max*deltaMS_max)*trapz(deltaBS,int2))^p;