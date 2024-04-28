clc;clear;close all;
%% Wire definitions
wireLength = 1;
a = 0.001;  %wire Radius
NoOfSegments= 20;
eps0 = 8.8542e-12;
delta = wireLength/NoOfSegments;
subPoints = 3;

basisPoints = 0:delta:wireLength;
a_sq = a*a;
%% Calculating Impedance matrix and Current along wire
Z = zeros;
for m = 1:NoOfSegments
    for n = 1:NoOfSegments
        l = sqrt(a_sq + (basisPoints(m)-basisPoints(n))^2);
        d_plus = l + delta/2;
        d_minus = l - delta/2;
        if m==n
            Z(m,n) = 2*log((delta/2 + sqrt(a_sq + (0.5*delta)^2))/a);
        elseif abs(m-n)<=2
            num = d_plus + sqrt(d_plus^2 + a_sq);
            denom = d_minus + sqrt(d_minus^2 + a_sq);
            Z(m,n) = log(num/denom);
        else
            Z(m,n) = log(d_plus/d_minus);
        end
    end
end

V = 4*pi*eps0*ones(NoOfSegments,1);

I = (Z\V).';
%% Plotting current along the wire
basisPoints_plot = linspace(0,1,NoOfSegments*50);
I_plot=zeros;
for m = 1:length(basisPoints_plot)
    for p = 1:length(basisPoints)-1
        if(basisPoints_plot(m)>=basisPoints(p) && basisPoints_plot(m)<=basisPoints(p+1))
            I_plot(m) = I(p);
        end
    end
end

plot(basisPoints_plot,I_plot*1e12,'k','LineWidth',1.4);
ylim([7 10.5]);
xlabel('Wire length (m)');
ylabel('Charge density (pC/m)');