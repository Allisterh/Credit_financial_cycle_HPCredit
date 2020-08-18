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
    A = arrayfun(@(x1,x2) stationaritytest(x1,x2), X, Y)
    B = double(A);
    % Generating the Z Data

    Z1 = X./(1 + abs(X))+Y./(1 + abs(Y));
    Z2 = -1.*(X./(1 + abs(X)).*(Y./(1 + abs(Y))));
    figure(1);              % Generating a new window to plot in.
    grid on
    hold on % Hold position should be at the start of plotting to overlap plots
    surf(X,Y,Z1,'FaceAlpha',0.5,'EdgeColor','none','FaceColor','red');             % The surface plotting function.
    surf(X,Y,Z2,'FaceAlpha',0.5,'EdgeColor','none','FaceColor',	'#4DBEEE');             % The surface plotting function.
    surf(X,Y,B,'FaceAlpha',0.5,'EdgeColor','none');             % The surface plotting function.
    
    zlabel('Transformed coefficient')
    xlabel('phi1-input')
    ylabel('phi2-input')
    legend('phi1','phi2','stationary')

    
%Morley
    
%Phi_y1    
    clc; clear all ;
    x = -5:0.1:5;             % The range of x values.
    y = -5:0.1:5;             % The range of y values.
    [X,Y] = meshgrid (x,y); % This generates the actual grid of x and y values.
    % Generating the Z Data
    aaa = X./(1 + abs(X));
    ccc = (1 - abs(aaa)).*Y./(1 + Y) + abs(aaa) - aaa^2;
    
    Z1 = 2*aaa;
    Z2 = -1*(aaa^2 + ccc); 
    
%     A = arrayfun(@(x1,x2) stationaritytest(x1,x2), X, Y)
    A2 = arrayfun(@(x1,x2) stationaritytest2(x1,x2), X, Y);
    B2 = double(A2);

%     figure(1);              % Generating a new window to plot in.
%     grid on
%     hold on % Hold position should be at the start of plotting to overlap plots
%     surf(X,Y,Z1,'FaceAlpha',0.5,'EdgeColor','none');             % The surface plotting function.
%     surf(X,Y,Z2,'FaceAlpha',0.5,'EdgeColor','none');             % The surface plotting function.
    
    figure(2);              % Generating a new window to plot in.
    grid on
    hold on % Hold position should be at the start of plotting to overlap plots
    surf(X,Y,Z1,'FaceAlpha',0.5,'EdgeColor','none','FaceColor','red');             % The surface plotting function.
    surf(X,Y,Z2,'FaceAlpha',0.5,'EdgeColor','none','FaceColor',	'#4DBEEE');             % The surface plotting function.
    surf(X,Y,B2,'FaceAlpha',0.5);             % The surface plotting function.

    zlabel('Transformed coefficient');
    xlabel('phi1-input');
    ylabel('phi2-input');
    legend('phi1','phi2','stationary') ;
    hold off
   
%=============%
%Inverse of a syms function
%=============%
syms y x;
f = exp(-x);
g=finverse(f);
g


%======================%
%Covariance matrix%
%======================%
clc; clear all ;

    ggg = [c11,0,0,0;
          0,c22,0,0;
          c31,0,c33,0;
          0,c42,0,c44];
        qqq = ggg*ggg.'
        
syms ggg c11 c22 c44 c55 c41 c52  c42   c21 c31 c32 c33;
    ggg = [c11,0,0,0,0,0;
          0,c22,0,0,0,0;
          0,0,0,0,0,0;
          c41,0,0,c44,0,0;
          0,c52,0,0,c55,0;
          0,0,0,0,0,0];
      
      ggg = [c11,0,0,0,0,0;
          0,c22,0,0,0,0;
          0,0,0,0,0,0;
          c41,0,0,c44,0,0;
          0,c52,0,0,c55,0;
          0,0,0,0,0,0];
      
      ggg
      ggg.'
    qqq = ggg*ggg.';
    qqq
    
      syms a11 a12 a21 a22;
      bb = [a11, a12;
            a21, a22];
      cc=bb.*bb';
      cc

    qqq = ggg*ggg';
    qqq
    
    
    ggg = [c11,0,0;
          c21,c22,0;
          c31,c32,c33];

    qqq = ggg*ggg.';
    qqq
    
    
    %=======================%
    %stationarity test
    %=======================%
    
Mdl = arima('AR',0.2,'MA',-0.5,'D',1,'Constant',0,...
'Variance',1.5);
T = 30;
rng(5);
Y = simulate(Mdl,T);
adftest(Y)

num = 1;
denom = LagOp([1 1.75 -0.76]);
quot = mrdivide(num,denom);

[r1,r2] = isStable(denom)
[r1,r2] = isStable(quot)

num = 1;
denom = LagOp([1 -0.9 ]);
quot = mrdivide(num,denom);

[r1,r2] = isStable(denom)

%https://www.mathworks.com/matlabcentral/answers/383108-use-meshgrid-to-evaluate-anonymous-function-which-accepts-vector-inputs

