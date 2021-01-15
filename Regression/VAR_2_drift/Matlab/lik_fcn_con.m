function val = lik_fcn_con(prmtr,y,T,START,prior)

%========================================================================%

  %Transform hyperparamters to impose constraints
  phi_y1 = prmtr(1);
  phi_y2 = prmtr(2);
  phi_yx = 0;
  
  phi_h1 = prmtr(3);
  phi_h2 = prmtr(4);
  phi_hx = 0;
  
  sig_nyy = prmtr(5)^2; % s.e. of HP permanent component
  sig_eyy = prmtr(6)^2; % s.e. of credit permanent component
  sig_nhh = prmtr(7)^2; % s.e. of the HP AR component
  sig_ehh = prmtr(8)^2; % s.e. of the credit AR component
  sig_eyeh = prmtr(9)*sqrt(sig_eyy*sig_ehh); % correlation of transitory components
  sig_nynh = prmtr(10)*sqrt(sig_nyy*sig_nhh); % correlation of permanent components

  mu = prmtr(11);
  
    F = [1,0,0,0,0,0; %Transition matrix
         0,phi_y1,phi_y2,0,phi_yx,0;
         0,1,0,0,0,0;
         0,0,0,1,0,0;
         0,phi_hx,0,0,phi_h1,phi_h2;
         0,0,0,0,1,0];

    Fstar = [phi_y1,phi_y2,phi_yx,0;
             1,0,0,0;
             phi_hx,0,phi_h1,phi_h2;
             0,0,1,0]; %Transition matrix of I(0) part, no trends;
    
    H = [1,1,0,0,0,0; %Measurement equation
        0,0,0,1,1,0];

    Q = [sig_nyy,0,0,sig_nynh,0,0; %Cov matrix
         0,sig_eyy,0,0,sig_eyeh,0;
         0, 0, 0, 0, 0, 0;
         sig_nynh, 0, 0, sig_nhh, 0, 0;
         0,sig_eyeh,0,0,sig_ehh,0;
         0, 0,0, 0, 0, 0];

    Qstar = [sig_eyy,0,sig_eyeh,0; %Cov matrix of I(0) part
             0,0,0,0;        
             sig_eyeh,0,sig_ehh,0;
             0,0,0,0];

    A = [0;0];

    beta_ll = [prior(1),prior(7),prior(7),prior(2),prior(8),prior(8)]'; %Starting values

    vecQstar = reshape(Qstar,[numel(Qstar),1]);
    vecP_ll = inv(eye(16) - kron(Fstar,Fstar))*vecQstar;
    
    %Var matrix of initial state vector
    P_ll = [prior(3),0,0,0,0,0;
            0,vecP_ll(1,1),0,0,vecP_ll(3,1),0;
            0,0,0,0,0,0;
            0,0,0,prior(4),0,0;
            0,vecP_ll(9,1),0,0,vecP_ll(11,1),0;
            0,0,0,0,0,0];
    
    muvec = [0,0,0,mu,0,0]'; %Drift vector
    
%========================================================================%
                
    lik_mat = zeros(T,1);
    
    %j_iter = 1;
    
    for j_iter = 1:T
        beta_tl  = muvec + F*beta_ll;
        P_tl  =  F*P_ll*F' + Q;

        vt = y(j_iter,1:2)' - H*beta_tl - A; %Prediction error

        ft =  H*P_tl*H'; %Variance of forecast error

        beta_tt = beta_tl + P_tl*H'*inv(ft)*vt;
        P_tt = P_tl - P_tl*H'*inv(ft)*H*P_tl;

        lik_mat(j_iter,1) = prior(5)*log(((2*pi)^2)*det(ft)) + prior(6)*vt'*inv(ft)*vt;

        beta_ll = beta_tt;
        P_ll = P_tt;
    end

    val = sum(lik_mat(START:T));

end