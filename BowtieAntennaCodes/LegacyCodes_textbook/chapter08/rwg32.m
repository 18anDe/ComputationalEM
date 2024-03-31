%RWG32 FREQUENCY LOOP
%   Calculates radiated field at a point as a function 
%   of frequency
%   The point is outside the metal surface
%   Uses the mesh file from RWG2, mesh2.mat, and
%   the file containing surface current coefficients,
%   current.mat, from RWG31 as inputs.
%
%   The following parameters need to be specified:
%   
%   Observation point           ObservationPoint[X; Y; Z] (m)
%
%   Copyright 2002 AEMM. Revision 2002/03/26
%   Chapter 8


clear all
%Load the data
load('mesh2');
load('current');
ObservationPoint=[0 0 1]'; %Observation point

%FRequency series
for FF=1:NumberOfSteps   
    FF
    omega=2*pi*f(FF);
    k=omega/c_;
    K=j*k;
    
    for m=1:EdgesTotal
        Point1=Center(:,TrianglePlus(m));
        Point2=Center(:,TriangleMinus(m));
        DipoleCenter(:,m)=0.5*(Point1+Point2);
        DipoleMoment(:,m)=EdgeLength(m)*CURRENT(m,FF)*(-Point1+Point2); 
    end    
    [EP,HP]=point(ObservationPoint,eta_,K,DipoleMoment,DipoleCenter);
    EField=sum(EP,2); HField=sum(HP,2);
    E(1:3,FF)=EField;
end

FileName=strcat('radiatedfield.mat'); 
save(FileName,'f','E','ObservationPoint');
plot(f,abs(E(2,:)),f,abs(E(1,:)),f,abs(E(3,:)))
grid on
xlabel('f,Hz')
ylabel('Electric field magnitude (three components, V/m)')
b=figure;
plot(f,unwrap(angle(E(2,:))))
grid on
xlabel('f,Hz')
ylabel('Phase of Ey, rad')