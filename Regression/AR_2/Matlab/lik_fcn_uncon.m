function val = lik_fcn_uncon(prmtr,y,T,START,prior)

%========================================================================%

  %Transform hyperparamters to impose constraints
  prmtr=trans_uncon(prmtr);
  phi_y1 = prmtr(1);
  phi_y2 = prmtr(2);
   
 
  sig_nyy = prmtr(3)^2; % s.e. of permanent component
  sig_eyy = prmtr(4)^2; % s.e. of  AR component

  
  mu = 0;
  
    F = [1,0,0; %Transition matrix
         0,phi_y1,phi_y2;
         0,1,0];

    Fstar = [phi_y1,phi_y2;
             1,0]; %Transition matrix of I(0) part, no trends;
    
    H = [1,1,0]; %Measurement equation

    Q = [sig_nyy,0,0; %Cov matrix
         0,sig_eyy,0;
         0, 0, 0];
 

    Qstar = [sig_eyy,0; %Cov matrix of I(0) part
             0,0];        

    A = 0;

    beta_ll = [prior(1),prior(5),prior(6)]'; %Starting values

    vecQstar = reshape(Qstar,[numel(Qstar),1]);
    vecP_ll = inv(eye(4) - kron(Fstar,Fstar))*vecQstar;
    
    %Var matrix of initial state vector
    P_ll = [prior(3),0,0;
            0,vecP_ll(1,1),0;
            0,0,0];
    
    muvec = [mu,0,0]'; %Drift vector
    
%========================================================================%
                
    lik_mat = zeros(T,1);
    
    %j_iter = 1;
    
    for j_iter = 1:T
        beta_tl  = muvec + F*beta_ll;
        P_tl  =  F*P_ll*F' + Q;

        vt = y(j_iter,1)' - H*beta_tl - A; %Prediction error

        ft =  H*P_tl*H'; %Variance of forecast error

        beta_tt = beta_tl + P_tl*H'*inv(ft)*vt;
        P_tt = P_tl - P_tl*H'*inv(ft)*H*P_tl;

    lik_mat(j_iter,1) = prior(3)*log(((2*pi)^2)*det(ft)) + prior(4)*vt'*inv(ft)*vt + ...
             prior(7)*(beta_tl(2)^2);   %GB

        beta_ll = beta_tt;
        P_ll = P_tt;
    end

    val = sum(lik_mat(START:T));

end