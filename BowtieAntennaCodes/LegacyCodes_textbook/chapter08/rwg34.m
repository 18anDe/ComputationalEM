%RWG34 Calculation of the received pulse using
%   antenna-to-antenna transfer function
%
%   Uses receivedfield.mat (created by RWG33)as an input
%
%   The following parameters need to be specified prior to 
%   calculations:
%   
%   Duration of the primary Gaussian voltage pulse (s)  Duration
%   Time offset to plot the received and the (delayed) 
%   transmitted pulse on the same figure including the 
%   "pulse tail"                                        Offset
%
%   Note: This script is not optimized and is rather slow
%
%   Copyright 2002 AEMM. Revision 2002/03/26 
%   Chapter 8

clear all
load('receivedfield')

%Identify the pulse parameters:

Duration=1.00e-9;   %in seconds
sigma=Duration/7;   %Gaussian pulse
Offset=500*sigma;   %Offset

%Sampling time
Ts=1/(2*max(f))     

%Length of DFT
N=2*NumberOfSteps;

%Obtain the pulse form
i=1; 
Time=-Duration; 
t(1)=Time;
p(1)=t(1)*exp(-t(1)^2/(2*sigma^2))/sigma^2;
while (Time<=Offset)
    i=i+1;
    t(i)=-Duration+i*Ts;
    Time=t(i);
    p(i)=t(i)*exp(-t(i)^2/(2*sigma^2))/sigma^2;
end
K=i;
p=p/max(p); %Normalization

%DFT
for n=1:N
    Int=0;
    for k=1:K
        Int=Int+p(k)*exp(-2*pi*j*n*k/N);
    end
    FDIRECT(n)=Int; %/pi
end

%To reproduce the pulse spectrum shown in 
%Fig.8.3b use the line:
%plot(f*sigma, abs(FDIRECT(1:length(f)))/pi); grid on;

%Spectrum continuation
for n=NumberOfSteps+1:N
    OutputVoltage(n)=conj(OutputVoltage(N+1-n));
    Impedance(n)    =conj(Impedance(N+1-n));
end

%Spectrum multiplication by transfer function
SPECTRUM=OutputVoltage.*FDIRECT;

%Alternatively use the line:
%SPECTRUM=OutputVoltage.*FDIRECT.*...
%(2*Impedance./(Impedance+50)).*(50./(Impedance+50));

%To check if the DFT/IDFT is correct use the line:
%SPECTRUM=FDIRECT;

%IDFT (received pulse)
for k=1:K
    Int=0;
    for n=1:N
        Int=Int+SPECTRUM(n)*exp(2*pi*j*n*k/N);
    end
    TRANSFORM(k)=real(Int)/N;
end

%Plot transmitted and received pulse
axis([-2 12 -1 1])
xlabel ('t, ns')
ylabel ('voltage, V')

hold on
p=p/max(abs(p));
plot(t*1e9,p);
%hold on;
plot(t*1e9,TRANSFORM,'.');
%hold on
plot(t*1e9,TRANSFORM);
grid on

%Second derivative of the Gaussian pulse -normalized
%(just in case...)
%t=[-5:0.1:5]; y=(-6*t+4*t.^3).*exp(-t.^2); plot(t,y)


