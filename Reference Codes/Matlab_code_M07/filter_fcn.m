function [beta_mat,fcst_mat] = filter_fcn(prmtr,y,T,START,prior)
    
    %Transform hyperparamters to impose constraints
    prmtr = trans(prmtr);
    
    phi_y1 = prmtr(1);
    phi_y2 = prmtr(2);
    phi_c1 = prmtr(3);
    phi_c2 = prmtr(4);
    mu = prmtr(5);
    sig_yy = prmtr(6)^2; %s.e. of income AR component
    sig_cc = prmtr(7)^2; %s.e. of consumption AR component
    sig_vv = prmtr(8)^2; %s.e. of the random walk component
    sig_yc = prmtr(9)*sqrt(sig_yy*sig_cc);
    sig_yv = prmtr(10)*sqrt(sig_yy*sig_vv);
    sig_cv = prmtr(11)*sqrt(sig_cc*sig_vv);
    cbar = prmtr(12);
    gamma = prmtr(13);

    F = [phi_y1,phi_y2,0,0,0; %Transition matrix
        1,0,0,0,0;
        0,0,phi_c1,phi_c2,0;
        0,0,1,0,0;
        0,0,0,0,1];

    Fstar = [phi_y1,phi_y2,0,0; %Transition matrix of I(0) part
            1,0,0,0;
            0,0,phi_c1,phi_c2;
            0,0,1,0];
    
    muvec = [0,0,0,0,mu]'; %Drift vector

    H = [1,0,0,0,1; %Measurement equation
        0,0,1,0,gamma];

    Q = [sig_yy,0,sig_yc,0,sig_yv; %Cov matrix
        0,0,0,0,0;
        sig_yc,0,sig_cc,0,sig_cv;
        0,0,0,0,0;
        sig_yv,0,sig_cv,0,sig_vv];

    Qstar = [sig_yy,0,sig_yc,0; %Cov matrix of I(0) part
            0,0,0,0;
            sig_yc,0,sig_cc,0;
            0,0,0,0];

    A = [0;cbar];

    beta_ll = [0,0,0,0,947.5]'; %Starting values

    vecQstar = reshape(Qstar,[numel(Qstar),1]);
    vecP_ll = inv(eye(16) - kron(Fstar,Fstar))*vecQstar;
    
    %Var matrix of initial state vector
    P_ll = [vecP_ll(1,1),0,vecP_ll(3,1),0,0;
            0,0,0,0,0;
            vecP_ll(9,1),0,vecP_ll(11,1),0,0;
            0,0,0,0,0;
            0,0,0,0,prior];
        
    beta_mat = zeros(T,5);
    fcst_mat = zeros(T,2);
    
    for j_iter = 1:T
        beta_tl  =  muvec + F*beta_ll;
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
        