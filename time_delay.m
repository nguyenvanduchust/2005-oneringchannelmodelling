function [tau_k] = time_delay(R, phi);

D = 10000*R;
c = 3*1e8;
tau_k = (R + sqrt(D^2 + R^2 + 2.*D.*R.*cos(phi))-D)/c;