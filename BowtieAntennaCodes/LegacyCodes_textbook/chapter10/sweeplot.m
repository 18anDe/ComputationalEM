%SWEEPLOT Plots the input impedance (resistance/reactance)
%   and reflection coefficient as functions of frequency
%
%   Uses current.mat as an input to load the impedance data
%   Run this script after RWG3
%
%   Copyright 2002 AEMM. Revision 2002/03/26 
%   Chapter 8/ Chapter 10


clear all
load('current.mat');

%plot impedance (real+imag)
a=figure
plot(f, real(Impedance),'.',f,imag(Impedance),'.');
hold on
plot(f, real(Impedance),f,imag(Impedance),'--');
xlabel ('Frequency, Hz')
ylabel('Input  resistance/reactance, Ohm')
title('Resistance-solid; reactance-dashed')
%axis([200e6 400e6 -100 250])
grid on

b=figure
%Plot input reflection coefficient
Gamma=(Impedance-50)./(Impedance+50);
Out=20*log10(abs(Gamma));
plot(f, Out);
xlabel ('Frequency, Hz')
ylabel ('Return loss, dB')
grid on
hold on

