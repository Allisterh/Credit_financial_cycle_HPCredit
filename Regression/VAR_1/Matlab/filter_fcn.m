function [beta_mat,fcst_mat] = filter_fcn(prmtr,y,T,START,prior)

%========================================================================%

  %Transform hyperparamters to impose constraints
     
  prmtr = trans(prmtr);
    
  phi_y1 = prmtr(1);
  phi_y2 = 0;
  phi_yx = 0;
  
  phi_h1 = prmtr(2);
  phi_h2 = 0;
  phi_hx = 0;
  
  sig_nyy = prmtr(3)^2; % s.e. of HP permanent component
  sig_eyy = prmtr(4)^2; % s.e. of credit permanent component
  sig_nhh = prmtr(5)^2; % s.e. of the HP AR component
  sig_ehh = prmtr(6)^2; % s.e. of the credit AR component
  sig_eyeh = prmtr(7)*sqrt(sig_eyy*sig_ehh);
  sig_nynh = prmtr(8)*sqrt(sig_nyy*sig_nhh);

  mu = 0;
  
    F = [1,0,0,0; %Transition matrix
         0,phi_y1,0,0;         
         0,0,1,0;
         0,0,0,phi_h1];

    Fstar = [phi_y1,0;
             0,phi_h1]; %Transition matrix of I(0) part, no trends;
    
    H = [1,1,0,0; %Measurement equation
        0,0,1,1];

    Q = [sig_nyy,0,sig_nynh,0; %Cov matrix
         0,sig_eyy,0,sig_eyeh;
         sig_nynh, 0,  sig_nhh, 0;
         0,sig_eyeh,sig_ehh,0];

    Qstar = [sig_eyy,sig_eyeh; %Cov matrix of I(0) part
             sig_eyeh,sig_ehh];

    A = [0;0];

    beta_ll = [prior(1),0,prior(2),0]'; %Starting values

    vecQstar = reshape(Qstar,[numel(Qstar),1]);
    vecP_ll = inv(eye(4) - kron(Fstar,Fstar))*vecQstar;
    
    %Var matrix of initial state vector
    P_ll = [prior(3),0,0,0;
            0,vecP_ll(1,1),0,vecP_ll(2,1);
            0,0,prior(4),0;
            0,vecP_ll(3,1),0,vecP_ll(4,1)];
    
    muvec = [0,0,mu,0]'; %Drift vector
    
%========================================================================%
               
    beta_mat = zeros(T,4);
    fcst_mat = zeros(T,2);
    
    for j_iter = 1:T
        beta_tl  = muvec + F*beta_ll;
        P_tl  =  F*P_ll*F' + Q;

        vt = y(j_iter,1:2)' - H*beta_tl - A; %Prediction error

        ft =  H*P_tl*H'; %Variance of forecast error

        beta_tt = beta_tl + P_tl*H'*inv(ft)*vt;
        P_tt = P_tl - P_tl*H'*inv(ft)*H*P_tl;

        beta_ll = beta_tt;
        P_ll = P_tt;
        
        fcst_mat(j_iter,:) = vt';
        beta_mat(j_iter,:) = beta_tt';
    end
        
end
        