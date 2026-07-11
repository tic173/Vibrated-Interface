function [Ac,v,zeta_all,w1,w2,w1p,w2p] = faradayFloquet_RT_boundary_GenEIG_cylindrical_pinned(a_guess,omega_star, R0, m, l, C, Bd, At, eta, n, varargin)
% s_star - growth rate
% omega_star - driving frequence
% beta_star - wave number
% n - order of expansion

mode_type       = 'SH';
g_sign = 1;


if nargin >=10
    mode_type = varargin{1};
end

if nargin >=11
    g_sign = varargin{2};
end


lb = 30;


%%


    if strcmpi(mode_type, 'SH')
        Ln  = 2*(n+1);
        B1   =  full( spdiags([ones(Ln,1) zeros(Ln,1) ones(Ln,1)],-1:1,Ln,Ln));

        % if n ==0
        %     B1 = eye(2);
        % end

        sj = -0.0+0.5*1i*omega_star;
     
    elseif strcmpi(mode_type, 'H')
        Ln  =  2*n+1;
        B1   =  full( spdiags([ones(Ln,1) zeros(Ln,1) ones(Ln,1)],-1:1,Ln,Ln));

        sj = 0+0i*omega_star;

    end

    beta_star = bessel_derivative_root(m, lb) ;

    Aj_P = zeros(Ln,lb);
    Aj_H = zeros(Ln,lb);

    for j = 1:lb

        [Aj_P(:,j),Aj_H(:,j)] = FD_RT_coeffs(n,sj,beta_star(j)/R0/1.0,omega_star,At,eta,C,Bd,mode_type,g_sign);
    end

    w = besselj(m,beta_star);


    for iter = 1:3

    lm_star1 = sqrt((1+At)/2/At/Bd/((g_sign+a_guess)));

    Im1 = besseli(m,R0/lm_star1);
    % Im_inv = besseli(m,lm_star/R0)
    Im_prime1  = 1/2*(besseli(m-1,R0/lm_star1)+besseli(m+1,R0/lm_star1));
    lambda1 = 2*lm_star1/R0./(1+beta_star.^2*lm_star1^2/R0^2)*Im_prime1/Im1.*(beta_star.^2./(beta_star.^2-m^2));

    % figure; plot(real(lambda1),'r-'); hold on; plot(imag(lambda1),'b-')

    % Sigma1 = 1*diag(lambda1./w)*ones(lb,1)*w;

     Sigma1 = 1*diag(lambda1)*ones(lb,1)*ones(size(w));

    I1 = eye(Ln);

    % Jm = kron(w,I1(1:end,:));

    Q = 1;

    A = Q'*blkdiag(diag(Aj_P(:))*(eye(lb*Ln)-kron(Sigma1,I1)*1)+diag(Aj_H(:)))*Q;

    B  = Q'*kron(eye(lb),B1)*Q;


    tol = 1e-10;
    [V,D]   =  eig(A,B);

    cond(A);
    D1       =  (diag(D));
    % 
    D       =  real(D1(abs(imag(D1))<tol*abs(real(D1))  & real(D1)>=0 ) );
    V       =  V(:,(abs(imag(D1))<tol*abs(real(D1))  & real(D1)>=0 ));


    [D,idx]       =  sort(D,'ascend','ComparisonMethod','real');

  
    Ac =  min(D(real(D)>=0));

    V = V(:,idx);

    v_weights = (1+reshape(ones(Ln,1).*lambda1,[],1));

      % plot(real(v_weights(n+2:2*(n+1):end)))

     % v = Q*V(:,1).*(1-reshape(ones(Ln,1).*lambda,[],1));
 

    v = (eye(lb*Ln)-kron(Sigma1,I1)*1)*Q*V(:,1);

    D(1)*9.81;

    % a_guess_all = real(sum(reshape(V(:,1),[],lb).*Aj_P,2))./real(sum(reshape(V(:,1),[],lb),2))/2;

    % a_guess_all = real(sum(reshape(V(:,1),[],lb).*(Aj_P*(1-Sigma)+Aj_H),2))./real(sum(reshape(V(:,1),[],lb),2))/2;

    % plot(-n-1:n,real(sum(reshape(V(:,1),Ln,lb).*(Aj_P*(1-Sigma)+Aj_H),2))./real(sum(reshape(V(:,1),Ln,lb),2))/2);

      % a_guess = a_guess_all(n+2)

       a_guess = -Ac;

      % if g_sign ==1
      %     a_guess = -Ac;
      % else
      % 
      %     a_guess = -Ac;
      %     % if g_sign+Ac<0
      %     %     a_guess = 0;
      %     % else
      %     %     a_guess = -1*max(Ac,1+1e-6);
      %     % end
      % 
      % end
         % a_guess = max(a_guess_all(a_guess_all<=g_sign));

    end

    %%

    zeta_all_p  = null((A-Ac*B));

    zeta_all = reshape(zeta_all_p.*v,Ln,[]);

    w1=0;w2=0;w1p=0;w2p=0;

    for j = 1:lb

        [w1_j,w2_j,w1p_j,w2p_j] = FD_RT_BCcoeffs(n,sj,beta_star(j),omega_star,At,eta,C,Bd,mode_type,g_sign,zeta_all(:,j));

        w1  = w1+w1_j;
        w2  = w2+w2_j;
        w1p = w1p+w1p_j;
        w2p = w2p+w2p_j;
    end

end

%%


function [An_P,An_H] = FD_RT_coeffs(n,s_star,beta_star,omega_star,At,eta,C,Bd,mode_type,g_sign)


if strcmpi(mode_type, 'SH')
    Ln  = -n-1:n;

elseif strcmpi(mode_type, 'H')
    Ln  =  -n:n;

end



xi = (1-At)/(1+At);

An_P = zeros(length(Ln),1);
An_H = zeros(length(Ln),1);

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

            if q1==1

            L6 = [1,0,1,0,0,0,0,0];

            else

            L6 = [(1+q1)*exp(-beta_star*2),(-1+q1),2*q1*exp(-1*beta_star*(q1+1)),0,0,0,0,0];

            end

            if real(q2)<1
                L7 = [0,0,0,0,1,exp(-beta_star*(2)),exp(-beta_star*(1-q2)),exp(-1*beta_star*(1+q2))];
            else
                L7 = [0,0,0,0,exp(beta_star*(1-q2)),exp(-beta_star*(1+q2)),1,exp(-2*beta_star*q2)];
            end

            if q2 ==1
                L8 = [0,0,0,0,0,1,0,1];

            else
                L8 = [0,0,0,0,(1-1/q2),(1+1/q2)*exp(-beta_star*2),0,2*exp(-1*beta_star*(1+q2))];
            end

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

    % if q1==1 && q2==1
    % 
    %               L1 = [1,1,0,0];
    %     L2 = [1,1,-1,-1];
    %     L3 = [1,-q1,1,q2];
    %     L4 = [2,(q1^2+1),-2*eta,-eta*(1+q2^2)];
    % 
    %     LHS = [ L1;L2;L3;L4];
    % 
    %     coeffs0 = null(LHS);
    %     coeffs = [coeffs0(1); 0; coeffs0(2);0; 0; coeffs0(3);0;coeffs0(4)];
    % 
    % end


    dwdz =  (coeffs(1)-coeffs(2)+coeffs(3)*q1-coeffs(4)*q1)*beta_star;

    d3wdz3 =  (coeffs(1)-coeffs(2)+coeffs(3)*q1^3-coeffs(4)*q1^3)*beta_star^3;

    An_P(n_idx) = 2*dwdz*((sn_star)/(beta_star^2)+3*C*(1-eta)*(1+At)/2/At)-2*C*(1-eta)*d3wdz3/(beta_star^2)*(1+At)/2/At;

    An_H(n_idx) = 2*(g_sign+beta_star^2*(1+At)/Bd/2/At);


    if nj ==0 && strcmpi(mode_type, 'H')
         An(n_idx) = (2./beta_star).*(g_sign*beta_star + (1+At)/2/Bd.*beta_star.^3);
    end

    n_idx = n_idx+1;
end


end




function [w1,w2,w1p,w2p] = FD_RT_BCcoeffs(n,s_star,beta_star,omega_star,At,eta,C,Bd,mode_type,g_sign,zeta_all)


% Nz = 10001;
%
% z = (linspace(-0.25,0.25,Nz))';

Nz = 20001;

alpha = 4.5;       % larger = more clustering near center

s = linspace(-1,1,Nz)';
z =  sinh(alpha*s) / sinh(alpha);


if strcmpi(mode_type, 'SH')
    Ln  = -n-1:n;

elseif strcmpi(mode_type, 'H')
    Ln  =  -n:n;

end

w1 = zeros(Nz,length(Ln));
w2 = zeros(Nz,length(Ln));
w1p = zeros(Nz,length(Ln));
w2p = zeros(Nz,length(Ln));

xi = (1-At)/(1+At);

n_idx = 1;

for nj = Ln

    q1 =  sqrt(1+(s_star+1i*nj*omega_star)/(C*beta_star^2));
    q2 =  sqrt(1+(s_star+1i*nj*omega_star)*xi/(eta*C*beta_star^2));

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

            RHS = [s_star+1i*nj*omega_star;0;0;0;0;0;0;0];

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

            RHS = [s_star+1i*nj*omega_star;0;0;0;0;0];

            coeffs0 = (LHS)\RHS;
            coeffs = [coeffs0(1); 0; coeffs0(2);0; coeffs0(3:end)];

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
        L3 = [1,q1,1,q2];
        L4 = [2,(q1^2+1),-2*eta,-eta*(1+q2^2)];

        LHS = [ L1;L2;L3;L4];

        RHS = [s_star+1i*nj*omega_star;0;0;0];

        coeffs0 = (LHS)\RHS;
        coeffs = [coeffs0(1); 0; coeffs0(2);0; 0; coeffs0(3);0;coeffs0(4)];

    end



    % coeffs_all(n_idx,:) = coeffs;

    w1(:,n_idx)   =  zeta_all(n_idx)*(coeffs(1)*exp(beta_star.*z) +coeffs(2)*exp(-beta_star.*z) +coeffs(3)*exp(beta_star*q1.*z) +coeffs(4)*exp(-beta_star*q1.*z));

    w1p(:,n_idx)   =   zeta_all(n_idx)*(beta_star*coeffs(1)*exp(beta_star.*z)  -beta_star*coeffs(2)*exp(-beta_star.*z) +beta_star*q1*coeffs(3)*exp(beta_star*q1.*z) -beta_star*q1*coeffs(4)*exp(-beta_star*q1.*z));


    coeffs(1)-coeffs(2)+q1*coeffs(3)-q1*coeffs(4)  -(coeffs(5)-coeffs(6)+q2*coeffs(7)-q2*coeffs(8)  );

    % coeffs(1)*exp(beta_star*-1)+coeffs(3)*exp(beta_star*q1*-1)
    %
    % coeffs(6)*exp(-beta_star*0.2442)+coeffs(8)*exp(-beta_star*q2*0.2442);


    if norm(coeffs(7))<1e-6
        w2(:,n_idx)   =  zeta_all(n_idx)*(coeffs(5)*exp(beta_star.*z) +coeffs(6)*exp(-beta_star.*z) +coeffs(8)*exp(-beta_star*q2.*z));

        w2p(:,n_idx)   =  zeta_all(n_idx)*( beta_star*coeffs(5)*exp(beta_star.*z)  -beta_star*coeffs(6)*exp(-beta_star.*z) -beta_star*q2*coeffs(8)*exp(-beta_star*q2.*z));

    else

        w2(:,n_idx)   =  zeta_all(n_idx)*(coeffs(5)*exp(beta_star.*z) +coeffs(6)*exp(-beta_star.*z) +coeffs(7)*exp(beta_star*q2.*z) +coeffs(8)*exp(-beta_star*q2.*z));

        w2p(:,n_idx)   =   zeta_all(n_idx)*(beta_star*coeffs(5)*exp(beta_star.*z)  -beta_star*coeffs(6)*exp(-beta_star.*z) +beta_star*q2*coeffs(7)*exp(beta_star*q2.*z) -beta_star*q2*coeffs(8)*exp(-beta_star*q2.*z));

    end


    % w2(:,n_idx)   =  coeffs(6)*exp(-beta_star.*z) +coeffs(8)*exp(-beta_star*q2.*z);

    % w2p(:,n_idx)   =  -beta_star*coeffs(6)*exp(-beta_star.*z) -beta_star*q2*coeffs(8)*exp(-beta_star*q2.*z);



    n_idx = n_idx+1;
end


end

