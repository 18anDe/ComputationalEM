%SWEEPLOT Plots the input impedance (resistance/reactance)
%   and reflection coefficient as functions of frequency
%
%   Uses current.mat as an input to load the impedance data
%   Run this script after RWG31
%
%   Copyright 2002 AEMM. Revision 2002/03/26 
%   Chapter 8

clear all
load('current.mat');

%Plot impedance (real+imag)
a=figure
plot(f, real(Impedance));
xlabel ('Frequency, Hz')
ylabel('Input  resistance, Ohm')
axis([0 5000e6 -500 500])
grid on
b=figure
plot(f, imag(Impedance));
xlabel ('Frequency, Hz')
ylabel('Input  reactance, Ohm')
axis([0 5000e6 -500 500])
grid on

%Plot reflection coefficient(return loss) versus 50 Ohm
c=figure
Gamma=(Impedance-50)./(Impedance+50);
Out=20*log10(abs(Gamma));
plot(f, Out);
xlabel ('Frequency, Hz')
ylabel('Reflection Coefficient, dB')
grid on

