function [s_star,zeta_all,w1,w2] = faradayFloquet_RT_GenEIG_cylindrical_pinned(Ac, omega_star, R0, m, l, C, Bd, At, eta, n, varargin)
% s_star - growth rate
% omega_star - driving frequence
% beta_star - wave number
% n - order of expansion

mode_type       = 'SH';
g_sign = 1;


if nargin >=11
    mode_type = varargin{1};
end

if nargin >=12
    g_sign = varargin{2};
end



lb = 30;




%%  Construct the matrix A

beta_star = bessel_derivative_root(m, lb) ;


if strcmpi(mode_type, 'SH')
    Ln  = 2*(n+1);
    B1  =  full( spdiags([ones(Ln,1) zeros(Ln,1) ones(Ln,1)],-1:1,Ln,Ln));

    B_sin = 1i*full( spdiags([-1*ones(Ln,1) zeros(Ln,1) 1*ones(Ln,1)],-1:1,Ln,Ln));

elseif strcmpi(mode_type, 'H')
    Ln  =  2*n+1;
    B1   =  full( spdiags([ones(Ln,1) zeros(Ln,1) ones(Ln,1)],-1:1,Ln,Ln));

end


[A0] = @(sj) FD_RT_coeffs(m,n,R0,sj,beta_star/R0,omega_star,At,eta,C,Bd,mode_type,g_sign);


B  = kron(eye(lb),B1);


%%    Solve 

fun = @(sj) det(((A0(sj))-Ac*B)/82);

opt=optimset('Maxiter',2000,'TolX',1e-10,'Tolfun',1e-10);
% 
% fun(1)
% 
% fun(1+1*beta_star(1)/2.5-1*beta_star(1)^2*C/2+sqrt(beta_star(1)+C*(beta_star(1)).^3)*1i)

s_star = fsolve(@(sj) fun(sj),1+1*beta_star(1)/2.5-1*beta_star(1)^2*C/2+sqrt(beta_star(1)+C*(beta_star(1)).^3)*1i, opt);

% det(((A0(1))-Ac*B*0)/85)
% 
% 
%  s_star = fsolve(@(sj) fun(sj),1, opt);


 % fun(s_star)

 % s_star = fsolve(@(sj) fun(sj), 1+1*beta_star/2.5-1*beta_star^2*C/1+sqrt(beta_star+C*(beta_star).^3)*1i*1.25, opt);

if imag(s_star)<1e-4 && real(s_star)<0
    s_star = fsolve(@(sj) fun(sj), abs(s_star), opt);
end


%%   Periodic Harmonic Components 

 zeta_all  = null(((A0(s_star))+Ac*B));
% 
 % [w1,w2] = FD_RT_BCcoeffs(n,s_star,beta_star,omega_star,At,eta,C,Bd,mode_type,g_sign);

 w1 = 0; w2 =0;


  % plot(t/(2*pi/omega_star),real(exp(growthrate_FD(j)*t).*transpose(sum(zeta_all(j,2:end).*exp(1i*omega_star*(-n:n).*t'),2)/sum(zeta_all(j,2:end).*exp(1i*omega_star*(-n:n).*(0+1*phase)),2))),'linewidth',2,'LineStyle','-','Color','r');


end



function [A] = FD_RT_coeffs(m,n,R0,s_star,beta_star_all,omega_star,At,eta,C,Bd,mode_type,g_sign)

xi = (1-At)/(1+At);


if strcmpi(mode_type, 'SH')
    Ln  = -n-1:n;

elseif strcmpi(mode_type, 'H')
    Ln  =  -n:n;

end

lb = 30;
lm_star = sqrt((1+At)/2/At/Bd/g_sign);

w = besselj(m,beta_star_all);
Im = besseli(m,R0/lm_star);
% Im_inv = besseli(m,lm_star/R0)
Im_prime  = 1/2*(besseli(m-1,R0/lm_star)+besseli(m+1,R0/lm_star));
lambda = 2*lm_star/R0./(1+beta_star_all.^2*lm_star^2/R0^2)*Im_prime/Im.*(beta_star_all.^2./(beta_star_all.^2-m^2))./w;

Sigma = 1*diag(lambda)*ones(lb,1)*w;


Aj_P = zeros(length(Ln),lb);
Aj_H = zeros(length(Ln),lb);


for j = 1:lb
    beta_star = beta_star_all(j);

    n_idx = 1;

    for nj = Ln

        sn_star = s_star+1i*nj*omega_star;

        q1 =  sqrt(1+(sn_star)/(C*beta_star^2));
        q2 =  sqrt(1+(sn_star)*xi/(eta*C*beta_star^2));


        if abs(real(exp(beta_star*(-1+q1))))>1e-15 && abs(real(exp(beta_star*(-1+q2))))>1e-15


            if abs(real(exp(beta_star*(-2))))>1e-20 && abs(real(exp(beta_star*(-1-q1))))>1e-20

                L1 = [1,1,1,1,0,0,0,0];
                L2 = [1,1,1,1,-1,-1,-1,-1];
                L3 = [1,-1,q1,-q1,-1,1,-q2,q2];
                L4 = [2,2,(q1^2+1),(q1^2+1),-2*eta,-2*eta,-eta*(1+q2^2),-eta*(1+q2^2)];

                if real(q1)<1
                    L5 = [exp(-beta_star*(2)),1,exp(-1*beta_star*q1-beta_star),exp(-beta_star+beta_star*q1),0,0,0,0];
                else
                    L5 = [exp(-beta_star*(1+q1)),exp(beta_star*(1-q1)),exp(-2*beta_star*q1),1,0,0,0,0];
                end

                L6 = [(1+q1)*exp(-beta_star*2),(-1+q1),2*q1*exp(-1*beta_star*(q1+1)),0,0,0,0,0];

                if real(q2)<1
                    L7 = [0,0,0,0,1,exp(-beta_star*(2)),exp(-beta_star*(1-q2)),exp(-1*beta_star*(1+q2))];
                else
                    L7 = [0,0,0,0,exp(beta_star*(1-q2)),exp(-beta_star*(1+q2)),1,exp(-2*beta_star*q2)];
                end

                L8 = [0,0,0,0,(1-1/q2),(1+1/q2)*exp(-beta_star*2),0,2*exp(-1*beta_star*(1+q2))];

                LHS = [ L1;L2;L3;L4;L5;L6;L7;L8];

                RHS = [sn_star;0;0;0;0;0;0;0];

                coeffs = (LHS)\RHS;

            else

                L1 = [1,1,0,0,0,0];
                L2 = [1,1,-1,-1,-1,-1];
                L3 = [1,q1,-1,1,-q2,q2];
                L4 = [2,(q1^2+1),-2*eta,-2*eta,-eta*(1+q2^2),-eta*(1+q2^2)];

                if real(q2)<1
                    L7 = [0,0,1,exp(-beta_star*(2)),exp(-beta_star*(1-q2)),exp(-1*beta_star*(1+q2))];
                else
                    L7 = [0,0,exp(beta_star*(1-q2)),exp(-beta_star*(1+q2)),1,exp(-2*beta_star*q2)];
                end

                L8 = [0,0,(1-1/q2),(1+1/q2)*exp(-beta_star*2),0,2*exp(-1*beta_star*(1+q2))];

                LHS = [ L1;L2;L3;L4;L7;L8];

                RHS = [sn_star;0;0;0;0;0];

                coeffs0 = (LHS)\RHS;
                coeffs = [coeffs0(1); 0; coeffs0(2);0; coeffs0(3:end)];

            end


            if  abs(real(exp(beta_star*(-1+q1))))>1e3 && abs(real(exp(beta_star*(-1+q2))))>1e3

                L1 = [1,1,1,0,0,0];
                L2 = [1,1,1,-1,-1,-1];
                L3 = [1,-1,q1,-1,1,q2];
                L4 = [0,0,1,0,0,-eta];

                % L5 = [exp(-beta_star*(2)),1,0,0,0,0];
                %
                % L6 = [(1+q1)*exp(-beta_star*2),(-1+q1),0,0,0,0];

                L5 = [exp(-beta_star*(2)),-exp(beta_star*(0)),q1*exp(-beta_star*q1-beta_star),0,0,0];

                L6 = [0,0,0, exp(beta_star*(0)),-exp(beta_star*(-2)),-q2*exp(-beta_star*q2-beta_star)];

                LHS = [ L1;L2;L3;L4;L5;L6];

                RHS = [sn_star;0;0;0;0;0];

                cond(LHS);

                coeffs0 = (LHS)\RHS;
                coeffs = [coeffs0(1); coeffs0(2); coeffs0(3);0; coeffs0(4); coeffs0(5);0;coeffs0(6)];


            end
            % eps = 1e-10;
            % [U,S,V] =  svd(LHS);
            % S = diag(S);
            %
            % idx = (S>eps);
            %
            % coeffs = (V(:,idx)*diag(1./S(idx))*U(:,idx)')*RHS;

        else
            % reduce to infinite height for large wave number

            L1 = [1,1,0,0];
            L2 = [1,1,-1,-1];
            L3 = [1,-q1,1,q2];
            L4 = [2,(q1^2+1),-2*eta,-eta*(1+q2^2)];

            LHS = [ L1;L2;L3;L4];

            RHS = [sn_star;0;0;0];

            coeffs0 = (LHS)\RHS;
            coeffs = [coeffs0(1); 0; coeffs0(2);0; 0; coeffs0(3);0;coeffs0(4)];

        end


        dwdz =  (coeffs(1)-coeffs(2)+coeffs(3)*q1-coeffs(4)*q1)*beta_star;

        d3wdz3 =  (coeffs(1)-coeffs(2)+coeffs(3)*q1^3-coeffs(4)*q1^3)*beta_star^3;

        Aj_P(n_idx,j) = 2*dwdz*((sn_star)/(beta_star^2)+3*C*(1-eta)*(1+At)/2/At)-2*C*(1-eta)*d3wdz3/(beta_star^2)*(1+At)/2/At;

        Aj_H(n_idx,j) = 2*(g_sign+beta_star^2*(1+At)/Bd/2/At);


        if nj ==0 && strcmpi(mode_type, 'H')
            An(n_idx) = (2./beta_star).*(g_sign*beta_star + (1+At)/2/Bd.*beta_star.^3);
        end

        n_idx = n_idx+1;
    end

    I1 = eye(length(Ln));

    A = blkdiag(diag(Aj_P(:))*(eye(lb*length(Ln))-kron(Sigma,I1)*1)+diag(Aj_H(:)));


end
end



function [w1,w2] = FD_RT_BCcoeffs(n,s_star,beta_star,omega_star,At,eta,C,Bd,mode_type,g_sign)


if strcmpi(mode_type, 'SH')
    Ln  = -n-1:n;

elseif strcmpi(mode_type, 'H')
    Ln  =  -n:n;

end

w1 = zeros(length(Ln),1);
w2 = zeros(length(Ln),1);

xi = (1-At)/(1+At);

n_idx = 1;

for nj = Ln

    q1 =  sqrt(1+(s_star+1i*nj*omega_star)/(C*beta_star^2));
    q2 =  sqrt(1+(s_star+1i*nj*omega_star)*xi/(eta*C*beta_star^2));


    L1 = [1,1,0,0];
    L2 = [1,1,-1,-1];
    L3 = [1,q1,1,q2];
    L4 = [2,(q1^2+1),-2*eta,-eta*(1+q2^2)];

    LHS = [ L1;L2;L3;L4];

    RHS = [s_star+1i*nj*omega_star;0;0;0];

    coeffs0 = (LHS)\RHS;
    coeffs = [coeffs0(1); 0; coeffs0(2);0; 0; coeffs0(3);0;coeffs0(4)];

   
    % coeffs_all(n_idx,:) = coeffs;

    % w1{n_idx}   =  @(z) coeffs(1)*exp(beta_star*z) +coeffs(2)*exp(-beta_star*z) +coeffs(3)*exp(beta_star*q1*z) +coeffs(4)*exp(-beta_star*q1*z);

        w1(n_idx)   =   (coeffs0(1) +coeffs0(2)) /s_star;

    % coeffs(1)*exp(beta_star*-1)+coeffs(3)*exp(beta_star*q1*-1) 
    % 
    % coeffs(6)*exp(beta_star*-1)+coeffs(8)*exp(beta_star*q2*-1) 
  
% w1{n_idx}(-1)
% 
% w1{n_idx}(0)
% 
% coeffs(5)
% 
% 
% coeffs(7)


    % w2{n_idx}   =  @(z) coeffs(5)*exp(beta_star*z) +coeffs(6)*exp(-beta_star*z) +coeffs(7)*exp(beta_star*q2*z) +coeffs(8)*exp(-beta_star*q2*z);

% w2{n_idx}(1)

    n_idx = n_idx+1;
end


end


