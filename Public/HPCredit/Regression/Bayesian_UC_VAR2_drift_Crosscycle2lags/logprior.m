function out=logprior(theta,F0,VF0,Q0,VQ0,rho0,Vrho0)


% % F~N(F0,VF0)
% F0=ones(2,1);
% VF0=eye(2)*0.2;
% %MU~N(MU0,VMU0);
% MU0=zeros(2,1);
% VMU0=eye(2);
% %1/R~Gamma(R0,VR0)
% R0=1;
% VR0=1;
% %1/Q(i,i)~Gamma(Q0,VQ0)
% Q0=0.1;
% VQ0=1;



out=0;
% % F~N(F0,VF0)
% F0=ones(2,1);
% VF0=eye(2)*0.2;
        a = F0;
        b = VF0;
        outi=log(mvnpdf(theta(1:8),a,b));
        out=out+outi;
        
% %MU~N(MU0,VMU0);
% MU0=zeros(2,1);
% VMU0=eye(2);
%         a = MU0;
%         b = VMU0;
%         outi=log(mvnpdf(theta(3),a,b));
%         out=out+outi;
%          
  %alpha has gamma prior with mean 3 and variance 1
  
% %1/R~Gamma(R0,VR0)
% R0=1;
% VR0=1;
%        mu=R0;
%         v=VR0;
%         b = v/mu;
%         a = mu/b;
%         outi=lpdfgam(theta(5),a,b);
%         out=out+outi;
%   
  % %1/Q(i,i)~Gamma(Q0,VQ0)
% Q0=0.1;
% VQ0=1;
        %mu=Q0;
        %v=VQ0;
        %b = v/mu;
        %a = mu/b;
        outi=(gampdf2(Q0(1),VQ0(1),1/(theta(9)^2*100))); 
        out=out+outi;
        outi=(gampdf2(Q0(2),VQ0(2),1/(theta(10)^2))); 
        out=out+outi;
        outi=(gampdf2(Q0(3),VQ0(3),1/(theta(11)^2*100))); 
        out=out+outi;
        outi=(gampdf2(Q0(4),VQ0(4),1/(theta(12)^2))); 
        out=out+outi;
        

        for i=13:14
          outi=log(normpdf(theta(i),rho0(i-12),Vrho0(i-12))); 
          out=out+outi;
        end

        %outi=gampdf1(VQ0,Q0,1/theta(3));
        %out=out+outi;
        
        %mu=Q0;
        %v=VQ0;
        %b = v/mu;
        %a = mu/b;
        %outi=gampdf1(VQ0,Q0,1/theta(4));
        %out=out+outi;
        
        
        
        
        
        
        
