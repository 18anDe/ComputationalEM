%VOLUMEMESH Creates solid volume metal meshes 
%   and their modificitations (see the text)
%
%   Particular cases:
%       - a cube
%       - junction of two rectangular waveguides
%       - mesh of open rectangular waveguides
%
%   The script is rather slow though
%
%   For some meshes and some observation angles the VIEWER function 
%   may display inaccurate results ("loose" existing triangles). 
%   Change the observation angle (rotate the mesh) in order to see 
%   the correct result
%   
%   Copyright 2002 AEMM. Revision 2002/04/23 
%   AppendixA    

clear all
warning off

%   The input data (see the text of Appendix A)
%   Only 1, 2, 4, 8, etc. elements per structure edge are allowed
Nx=4; Ny=4;
Size=0.5;
a=1; b=1;
h=0.5; H=1.0;


ElemArray   =[];
PointArray  =[];
Shift       =0;
N=Nx;
M=Ny;

% a sides for both Z
for k=0:1
    for i=1:N
   	    for j=1:M   
   		    X1 =(i-1)*(a+h); 
   		    X2 = i*a+(i-1)*h;
   		    Y1 = b*j+h*(j-1);
  		    Y2 = b*j+h*j;
   		    Z1 =k*H;   
   		    Z2 =k*H;
   		    [p, t] =strip_(X1, X2, Y1, Y2, Z1, Z2, Size);
            PointArray=[PointArray p];
            ElemArray=[ElemArray t+Shift];
      	    Shift=Shift+length(p);
        end
    end
end

% b sides for both Z
for k=0:1
    for i=1:N-1
        for j=1:M+1   
   		    X1 = i*a+(i-1)*h; 
   		    X2 = i*a+i*h;
   		    Y1 = b*(j-1)+h*(j-1);
   		    Y2 = b*j+h*(j-1);
   		    Z1 = k*H;   
   		    Z2 = k*H;
   		    [p, t] =strip_(X1, X2, Y1, Y2, Z1, Z2, Size);
            PointArray=[PointArray p];
            ElemArray=[ElemArray t+Shift];
      	    Shift=Shift+length(p);      
        end
    end
end

% h sides for both Z
for k=0:1
    for i=1:N-1
   	    for j=1:M   
   		    X1 = i*a+(i-1)*h; 
   		    X2 = i*a+i*h;
   		    Y1 = b*j+h*(j-1);
   		    Y2 = b*j+h*j;
   		    Z1 = k*H;   
   		    Z2 = k*H;
   		    [p, t] =strip_(X1, X2, Y1, Y2, Z1, Z2, Size);
            PointArray=[PointArray p];
            ElemArray=[ElemArray t+Shift];
      	    Shift=Shift+length(p);      
        end
    end
end

% Vertical sides a
for i=1:N
    for j=1:M   
   	    X1 = (i-1)*a+(i-1)*h; 
   	    X2 = i*a+(i-1)*h;
   	    Y1 = b*j+h*(j-1);
   	    Y2 = Y1;
   	    Z1 = 0;   
   	    Z2 = H;
        [p, t] =strip_(X1, X2, Y1, Y2, Z1, Z2, Size);
        PointArray=[PointArray p];
        ElemArray=[ElemArray t+Shift];
        Shift=Shift+length(p);      
    end
end

% Vertical sides a+h
for i=1:N
    for j=1:M   
   	    X1 = (i-1)*a+(i-1)*h; 
   	    X2 = i*a+(i-1)*h;
   	    Y1 = b*j+h*j;
   	    Y2 = Y1;
   	    Z1 = 0;   
   	    Z2 = H;
        [p, t] =strip_(X1, X2, Y1, Y2, Z1, Z2, Size);
        PointArray=[PointArray p];
        ElemArray=[ElemArray t+Shift];
        Shift=Shift+length(p);      
    end
end

% Vertical sides b
for i=1:N-1
    for j=1:M+1   
   	    X1 = i*a+(i-1)*h; 
   	    X2 = X1;
   	    Y1 = b*(j-1)+h*(j-1);
   	    Y2 = b*j+h*(j-1);
   	    Z1 = 0;   
   	    Z2 = H;
   	    [p, t] =strip_(X1, X2, Y1, Y2, Z1, Z2, Size);
        PointArray=[PointArray p];
        ElemArray=[ElemArray t+Shift];
        Shift=Shift+length(p);      
    end
end

% Vertical sides b+h
for i=1:N-1
    for j=1:M+1   
   	    X1 = i*a+i*h; 
   	    X2 = X1;
   	    Y1 = b*(j-1)+h*(j-1);
   	    Y2 = b*j+h*(j-1);
   	    Z1 = 0;   
   	    Z2 = H;
        [p, t] =strip_(X1, X2, Y1, Y2, Z1, Z2, Size);
        PointArray=[PointArray p];
        ElemArray=[ElemArray t+Shift];
        Shift=Shift+length(p);      
    end
end

%Eliminate this block to obtain the mesh of open waveguides
% H top
for i=1:N-1
    X1 = i*a+(i-1)*h; 
   	X2 = i*a+i*h;
   	Y1 = (M+1)*b+M*h;
   	Y2 = Y1;
   	Z1 = 0;   
   	Z2 = H;
   	[p, t] =strip_(X1, X2, Y1, Y2, Z1, Z2, Size);
    PointArray=[PointArray p];
    ElemArray=[ElemArray t+Shift];
    Shift=Shift+length(p);      
end

% H bottom
for i=1:N-1
    X1 = i*a+(i-1)*h; 
   	X2 = i*a+i*h;
   	Y1 = 0;
   	Y2 = Y1;
   	Z1 = 0;   
   	Z2 = H;
   	[p, t] =strip_(X1, X2, Y1, Y2, Z1, Z2, Size);
    PointArray=[PointArray p];
    ElemArray=[ElemArray t+Shift];
    Shift=Shift+length(p);      
end

% H left
for j=1:M
    X1 = 0; 
   	X2 = X1;
   	Y1 = j*b+(j-1)*h;
   	Y2 = j*b+j*h;
   	Z1 = 0;   
   	Z2 = H;
   	[p, t] =strip_(X1, X2, Y1, Y2, Z1, Z2, Size);
    PointArray=[PointArray p];
    ElemArray=[ElemArray t+Shift];
    Shift=Shift+length(p);      
end

% H right
for j=1:M
    X1 = N*a+(N-1)*h; 
   	X2 = X1;
   	Y1 = j*b+(j-1)*h;
   	Y2 = j*b+j*h;
   	Z1 = 0;   
   	Z2 = H;
   	[p, t] =strip_(X1, X2, Y1, Y2, Z1, Z2, Size);
    PointArray=[PointArray p];
    ElemArray=[ElemArray t+Shift];
    Shift=Shift+length(p);      
end
%Eliminate this block to obtain the mesh of open waveguides


p               =PointArray;
t               =ElemArray;
t(4,:)          =1;
PointNumber     =length(PointArray);
TrianglesTotal  =length(ElemArray);

MeanX=mean(PointArray(1,:));
MeanY=mean(PointArray(2,:));
MeanZ=mean(PointArray(3,:));
PointArray(1,:)=PointArray(1,:)-MeanX;
PointArray(2,:)=PointArray(2,:)-MeanY;
PointArray(3,:)=PointArray(3,:)-MeanZ;

% Reduce number of points
epsilon=1e-6;
P(1:3,1)=p(1:3,1);
INDEX=1;
for L=2:PointNumber
   Indicator=0;
   Point=p(1:3,L);   
   for k=1:INDEX
       if(norm(Point-P(1:3,k))<epsilon)
           Indicator=1;
       end
   end
   if (Indicator==0)
       P=[P Point];
       INDEX=INDEX+1;
   end      
end

% Re-enumerate triangles
for L=1:TrianglesTotal;
    n1=t(1,L);n2=t(2,L); n3=t(3,L);
    for k=1:length(P)
        if(norm(P(1:3,k)-p(1:3,n1))<epsilon)
            Number1=k;
        end
        if(norm(P(1:3,k)-p(1:3,n2))<epsilon)
            Number2=k;
        end
        if(norm(P(1:3,k)-p(1:3,n3))<epsilon)
            Number3=k;
        end
    end
    T(1:3,L)=[Number1; Number2; Number3];
    T(4,  L)=1;
end

clear t, p;
p=P; t=T;

% Save result
save mesh  p t;
viewer mesh