function [Ac] = faradayFloquet_RT_boundary_GenEIG( omega_star, k_star, C, Bd, At, eta, n, varargin)
% s_star - growth rate
% omega_star - driving frequence
% k_star - wave number
% n - order of expansion

mode_type       = 'SH';
g_sign = 1;


if nargin >=8
    mode_type = varargin{1};
end

if nargin >=9
    g_sign = varargin{2};
end


%%

    if strcmpi(mode_type, 'SH')
        Ln  = 2*(n+1);
        B   =  full( spdiags([ones(Ln,1) zeros(Ln,1) ones(Ln,1)],-1:1,Ln,Ln));

        sj = 0.5*1i*omega_star;
     
    elseif strcmpi(mode_type, 'H')
        Ln  =  2*n+1;
        B   =  full( spdiags([ones(Ln,1) zeros(Ln,1) ones(Ln,1)],-1:1,Ln,Ln));

        sj = 0;

    end

         A0 = FD_RT_coeffs(n,sj,k_star,omega_star,At,eta,C,Bd,mode_type,g_sign);


         tol = 1e-6;
    [V,D]   =  eig(diag(A0),B);
    D1       =  (diag(D));

     D       =  real(D1(abs(imag(D1))<tol*abs(real(D1))));

    D       =  sort(D,'ascend','ComparisonMethod','real');

     Ac =  min(D(real(D)>0));



end


function [An] = FD_RT_coeffs(n,s_star,k_star,omega_star,At,eta,C,Bd,mode_type,g_sign)


if strcmpi(mode_type, 'SH')
    Ln  = -n-1:n;

elseif strcmpi(mode_type, 'H')
    Ln  =  -n:n;

end


xi = (1-At)/(1+At);

An = zeros(length(Ln),1);

n_idx = 1;

for nj = Ln

    q1 =  sqrt(1+(s_star+1i*nj*omega_star)/(C*k_star^2));
    q2 =  sqrt(1+(s_star+1i*nj*omega_star)*xi/(eta*C*k_star^2));

     
    if abs(real(exp(k_star*(-1+q1))))>1e-15 && abs(real(exp(k_star*(-1+q2))))>1e-15


        if abs(real(exp(k_star*(-2))))>1e-20 && abs(real(exp(k_star*(-1-q1))))>1e-20

            L1 = [1,1,1,1,0,0,0,0];
            L2 = [1,1,1,1,-1,-1,-1,-1];
            L3 = [1,-1,q1,-q1,-1,1,-q2,q2];
            L4 = [2,2,(q1^2+1),(q1^2+1),-2*eta,-2*eta,-eta*(1+q2^2),-eta*(1+q2^2)];

            if real(q1)<1
                L5 = [exp(-k_star*(2)),1,exp(-1*k_star*q1-k_star),exp(-k_star+k_star*q1),0,0,0,0];
            else
                L5 = [exp(-k_star*(1+q1)),exp(k_star*(1-q1)),exp(-2*k_star*q1),1,0,0,0,0];
            end

            L6 = [(1+q1)*exp(-k_star*2),(-1+q1),2*q1*exp(-1*k_star*(q1+1)),0,0,0,0,0];

            if real(q2)<1
                L7 = [0,0,0,0,1,exp(-k_star*(2)),exp(-k_star*(1-q2)),exp(-1*k_star*(1+q2))];
            else
                L7 = [0,0,0,0,exp(k_star*(1-q2)),exp(-k_star*(1+q2)),1,exp(-2*k_star*q2)];
            end

            L8 = [0,0,0,0,(1-1/q2),(1+1/q2)*exp(-k_star*2),0,2*exp(-1*k_star*(1+q2))];

            LHS = [ L1;L2;L3;L4;L5;L6;L7;L8];

            RHS = [s_star+1i*nj*omega_star;0;0;0;0;0;0;0];

            coeffs = (LHS)\RHS;

        else

            L1 = [1,1,0,0,0,0];
            L2 = [1,1,-1,-1,-1,-1];
            L3 = [1,q1,-1,1,-q2,q2];
            L4 = [2,(q1^2+1),-2*eta,-2*eta,-eta*(1+q2^2),-eta*(1+q2^2)];

            if real(q2)<1
                L7 = [0,0,1,exp(-k_star*(2)),exp(-k_star*(1-q2)),exp(-1*k_star*(1+q2))];
            else
                L7 = [0,0,exp(k_star*(1-q2)),exp(-k_star*(1+q2)),1,exp(-2*k_star*q2)];
            end

            L8 = [0,0,(1-1/q2),(1+1/q2)*exp(-k_star*2),0,2*exp(-1*k_star*(1+q2))];

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
        L3 = [1,-q1,1,q2];
        L4 = [2,(q1^2+1),-2*eta,-eta*(1+q2^2)];

        LHS = [ L1;L2;L3;L4];

        RHS = [s_star+1i*nj*omega_star;0;0;0];

        coeffs0 = (LHS)\RHS;
        coeffs = [coeffs0(1); 0; coeffs0(2);0; 0; coeffs0(3);0;coeffs0(4)];

    end


    dwdz =  (coeffs(1)-coeffs(2)+coeffs(3)*q1-coeffs(4)*q1)*k_star;

    d3wdz3 =  (coeffs(1)-coeffs(2)+coeffs(3)*q1^3-coeffs(4)*q1^3)*k_star^3;

    An(n_idx) = 2*dwdz*((s_star+1i*nj*omega_star)/(k_star^2)+3*C*(1-eta)*(1+At)/2/At)-2*C*d3wdz3/(k_star^2)*(1+At)/2/At+2*(g_sign+k_star^2*(1+At)/Bd/2);

    if nj ==0 && strcmpi(mode_type, 'H')
         An(n_idx) = (2./k_star).*(g_sign*k_star + (1+At)/2/Bd.*k_star.^3);
    end

    n_idx = n_idx+1;
end


end


