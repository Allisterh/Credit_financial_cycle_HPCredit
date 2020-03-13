function val = lik_fcn(prmtr,y,T,START,t_h_prior, t_c_prior)

    %Transform hyperparamters to impose constraints
  prmtr = trans(prmtr);
    
  phi_h11 = prmtr(1);
  phi_h12 = prmtr(2);
  phi_h21 = prmtr(3);
  phi_h22 = prmtr(4);
  
  phi_c11 = prmtr(5);
  phi_c12 = prmtr(6);
  phi_c21 = prmtr(7);
  phi_c22 = prmtr(8);
  
  mu_h = prmtr(9);
  mu_c = prmtr(10);
    
  sig_nhh = prmtr(11)^2; % s.e. of HP permanent component
  sig_ncc = prmtr(12)^2; % s.e. of credit permanent component
  sig_ehh = prmtr(13)^2; % s.e. of the HP AR component
  sig_ecc = prmtr(14)^2; % s.e. of the credit AR component
  sig_nhnc = prmtr(15)*sqrt(sig_nhh*sig_ncc);
  sig_ehec = prmtr(16)*sqrt(sig_ehh*sig_ecc);

  

    F = [1,0,0,0,0,0; %Transition matrix
         0,phi_h11,phi_h12,0,phi_h21,phi_h22;
         0,1,0,0,0,0;
         0,0,0,1,0,0;
         0,phi_c11,phi_c12,0,phi_c21, phi_c22;
         0,0,0,0,1,0];

    Fstar = [phi_h11,phi_h12,phi_h21,phi_h22;
            1,0,0,0;
            phi_c11,phi_c12,phi_c21, phi_c22;
            0,0,1,0]; %Transition matrix of I(0) part];
    
    muvec = [mu_h,0,0,mu_c,0,0]'; %Drift vector

    H = [1,1,0,0,0,0; %Measurement equation
        0,0,0,1,1,0];

    Q = [sig_nhh,0,0,sig_nhnc,0,0; %Cov matrix
        0, sig_ehh, 0, 0, sig_ehec, 0;
        0,0,0,0,0,0;
        sig_nhnc, 0, 0, sig_ncc, 0, 0;
        0, sig_ehec, 0, 0, sig_ecc, 0;
        0,0,0,0,0,0];

    Qstar = [sig_ehh,0,sig_ehec,0; %Cov matrix of I(0) part
            0,0,0,0;
            sig_ehec, 0,sig_ecc, 0;
            0,0,0,0];

    A = [0;0];

    beta_ll = [t_h_prior,0,0,t_c_prior,0,0]'; %Starting values

    vecQstar = reshape(Qstar,[numel(Qstar),1]);
    vecP_ll = inv(eye(16) - kron(Fstar,Fstar))*vecQstar;
    
    %Var matrix of initial state vector
    P_ll = [100,0,0,50,0,0;
            0,vecP_ll(1,1),0,0,vecP_ll(3,1),0;
            0,0,0,0,0,0;
            70,0,0,200,0,0;
            0,vecP_ll(9,1),0,0,vecP_ll(11,1),0;
            0,0,0,0,0,0];

    lik_mat = zeros(T,1);
    
    for j_iter = 1:T
        beta_tl  =  muvec + F*beta_ll;
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