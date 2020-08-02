function val = lik_fcn(prmtr,y,T,START,prior)

  %Transform hyperparamters to impose constraints
     
  prmtr = trans(prmtr);
    
  phi_y1 = prmtr(1);
  phi_yh = prmtr(2);
  phi_h1 = prmtr(3);
  phi_hy = prmtr(4);
  
  sig_nyy = prmtr(5)^2; % s.e. of HP permanent component
  sig_nhh = prmtr(6)^2; % s.e. of credit permanent component
  sig_eyy = prmtr(7)^2; % s.e. of the HP AR component
  sig_ehh = prmtr(8)^2; % s.e. of the credit AR component
  sig_nynh = prmtr(9)*sqrt(sig_nyy*sig_nhh);
  sig_eyeh = prmtr(10)*sqrt(sig_eyy*sig_ehh);
  
  sig_wyy = prmtr(11);
  sig_whh = prmtr(12);

    F = [1,1,0,0,0,0; %Transition matrix
         0,1,0,0,0,0;
         0,0,phi_y1,0,0,phi_yh;
         0,0,0,1,1,0;
         0,0,0,0,1,0;
         0,0,phi_hy,0,0,phi_h1];

    Fstar = [phi_y1,phi_yh;
            phi_hy,phi_h1]; %Transition matrix of I(0) part];
    
    H = [1,0,1,0,0,0; %Measurement equation
        0,0,0,1,0,1];

    Q = [sig_nyy,0,0,sig_nynh,0,0; %Cov matrix
        0,sig_wyy,0,0,0,0;
        0, 0,sig_eyy,0, 0, sig_eyeh;
        sig_nynh, 0, 0, sig_nhh, 0, 0;
        0,0,0,0,sig_whh,0;
        0, 0,sig_eyeh, 0, 0, sig_ehh];

    Qstar = [sig_eyy,sig_eyeh; %Cov matrix of I(0) part
            sig_eyeh,sig_ehh];

    A = [0;0];

    beta_ll = [prior(1),0,0,prior(2),0,0]'; %Starting values

    vecQstar = reshape(Qstar,[numel(Qstar),1]);
    vecP_ll = inv(eye(4) - kron(Fstar,Fstar))*vecQstar;
    
    %Var matrix of initial state vector
    P_ll = [prior(3),0,0,prior(7),0,0;
            0,prior(4),0,0,0,0;
            0,0,vecP_ll(1,1),0,0,vecP_ll(2,1);
            prior(7),0,0,prior(5),0,0;
            0,0,0,0,prior(6),0;
            0,0,vecP_ll(3,1),0,0,vecP_ll(4,1)];
                
    lik_mat = zeros(T,1);
    
    %j_iter = 1;
    
    for j_iter = 1:T
        beta_tl  =  F*beta_ll;
        P_tl  =  F*P_ll*F' + Q;

        vt = y(j_iter,1:2)' - H*beta_tl - A; %Prediction error

        ft =  H*P_tl*H'; %Variance of forecast error

        beta_tt = beta_tl + P_tl*H'*inv(ft)*vt;
        P_tt = P_tl - P_tl*H'*inv(ft)*H*P_tl;

        lik_mat(j_iter,1) = 0.5*log(((2*pi)^2)*det(ft)) + 0.5*vt'*inv(ft)*vt;

        beta_ll = beta_tt;
        P_ll = P_tt;
    end

    val = sum(lik_mat(START:T));

end