%SWEEPLOT Plots the input impedance - 
%   resistance and reactance - as functions of frequency
%
%   Uses current.mat as an input to load the impedance data
%
%   Copyright 2002 AEMM. Revision 2002/03/26 Chapter 9

clear all
load current

%plot impedance (real+imag)
a=figure
plot(f, real(Impedance),f,imag(Impedance),'--');
xlabel ('Frequency, Hz')
ylabel('Input  resistance/reactance, Ohm')
title('Resistance-solid; reactance-dashed')
axis([0 5e8 -1000 1500])
grid on
