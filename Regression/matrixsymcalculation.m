syms a b c d a11 a12 a13 a14 a21 a22 a23 a24

aa = [a,0,0,b;
     0,0,0,0;
     0,0,0,0;
     c,0,0,d];

bb= aa*aa.';
bb

aa = [a11,a12,a13,0;
      0,1,0,0;
      a21,0,a23,a24;
      0,0,0,1];
  
bb = kron(aa,aa.');
bb
size(bb)


    

    

    
    %=============%
    % 3D Graph
    %=============%
    
    clc; clear all ;
    x = -20:0.25:20;             % The range of x values.
    y = -20:0.25:20;             % The range of y values.
    [X,Y] = meshgrid (x,y); % This generates the actual grid of x and y values.
    % Generating the Z Data

    Z = X./(1 + abs(X))+Y./(1 + abs(Y));
%     Z2 = -1* (aaa*ccc);
    figure(1);              % Generating a new window to plot in.
    surf(X,Y,Z)             % The surface plotting function.
    
%Phi_y2    
    clc; clear all ;
    x = -20:0.25:20;             % The range of x values.
    y = -20:0.25:20;             % The range of y values.
    [X,Y] = meshgrid (x,y); % This generates the actual grid of x and y values.
    % Generating the Z Data
    
    Z = -1.*(X./(1 + abs(X)).*(Y./(1 + abs(Y))));
%     Z2 = -1* (aaa*ccc);
    figure(1);              % Generating a new window to plot in.
    surf(X,Y,Z)             % The surface plotting function.
    
%Morley
    
%Phi_y1    
    clc; clear all ;
    x = -20:0.25:20;             % The range of x values.
    y = -20:0.25:20;             % The range of y values.
    [X,Y] = meshgrid (x,y); % This generates the actual grid of x and y values.
    % Generating the Z Data
    aaa = 2.*X./(1 + abs(X));
    ccc = (1 - abs(aaa)).*Y./(1 + Y) + abs(aaa) - aaa^2;
    
    Z1 = 2*aaa;
    Z2 = -1*(aaa^2 + ccc); 
%     Z2 = -1* (aaa*ccc);
    figure(1);              % Generating a new window to plot in.
    grid on
    hold on % Hold position should be at the start of plotting to overlap plots
    surf(X,Y,Z1);             % The surface plotting function.
    surf(X,Y,Z2,'FaceAlpha',0.5,'EdgeColor','none');             % The surface plotting function.
    

    
x= linspace(-20,20,100);
y = 2*x./(1 + abs(x));
plot(x,y)
y1 = 2*1.5./(1 + abs(1.5))
aaa = 1.5/2

x2=linspace(-20,20,100);
c = (1 - abs(aaa))*x2./(1 + abs(x2)) + abs(aaa) - aaa^2;
y2 = -1* (aaa^2 + c);
plot(x2,y2)

clear, clc
x2=-20;
aaa = 1.5/2;
c2 = (1 - abs(aaa))*x2./(1 + abs(x2)) + abs(aaa) - aaa^2;
y3= -1* (aaa^2 + c2)

c0= linspace(-20,20,n);
aaa = c0(1)./(1 + abs(c0(1)));
ccc = (1 - abs(aaa))*c0(2)./(1 + abs(c0(2))) + abs(aaa) - aaa^2;

c1(1) = 2*aaa;
c1(2) = -1* (aaa^2 + ccc);
    