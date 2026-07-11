function [An] = FD_RT_coeffs(n,s_star,k_star,omega_star,At,eta,C,Bd,mode_type)


if strcmpi(mode_type, 'SH')
    Ln  = -n-1:n;
    % B   =  full( spdiags([ones(Ln,1) zeros(Ln,1) ones(Ln,1)],-1:1,Ln,Ln));

elseif strcmpi(mode_type, 'H')
    Ln  =  -n:n;
    % B   =  full( spdiags([ones(Ln,1) zeros(Ln,1) ones(Ln,1)],-1:1,Ln,Ln));

end


xi = (1-At)/(1+At);

An = zeros(length(Ln),1);

n_idx = 1;

for nj = Ln

    q1 =  sqrt(1+(s_star+1i*nj*omega_star)/(C*k_star^2));
    q2 =  sqrt(1+(s_star+1i*nj*omega_star)*xi/(eta*C*k_star^2));

    LHS = [ 1,1,1,1,0,0,0,0;...
        1,1,1,1,-1,-1,-1,-1;...
        1,-1,q1,-q1,-1,1,-q2,q2;...
        2,2,(q1^2+1),(q1^2+1),-2*eta,-2*eta,-eta*(q2^2+1),-eta*(q2^2+1);...
        exp(-k_star*(1+q1)),exp(k_star*(1-q1)),exp(-2*k_star*q1),1,0,0,0,0;...
        (1+q1)*exp(-k_star*(1)),(-1+q1)*exp(k_star*(1)),2*q1*exp(-1*k_star*q1),0,0,0,0,0;...
        0,0,0,0,exp(k_star*(1-q2)),exp(-k_star*(1+q2)),1,exp(-2*k_star*q2);...
        0,0,0,0,(q2-1)*exp(k_star*(1)),(q2+1)*exp(-k_star*(1)),0,2*q2*exp(-1*k_star*q2)...
        ];

    RHS = [s_star+1i*nj*omega_star;0;0;0;0;0;0;0];


    eps = 1e-12;
    [U,S,V] =  svd(LHS);
    S = diag(S);

    idx = (S>eps);

    coeffs = (V(:,idx)*diag(1./S(idx))*U(:,idx)')*RHS;

    % coeffs = (LHS)\RHS;

    dwdz =  (coeffs(1)-coeffs(2)+coeffs(3)*q1-coeffs(4)*q1)*k_star;

    d3wdz3 =  (coeffs(1)-coeffs(2)+coeffs(3)*q1^3-coeffs(4)*q1^3)*k_star^3;

    An(n_idx) = 2*dwdz*((s_star+1i*nj*omega_star)/(k_star^2)+3*C*(1-eta)*(1+At)/2/At)-2*C*d3wdz3/(k_star^2)*(1+At)/2/At+2*(1+k_star^2*(1+At)/Bd/2);

    n_idx = n_idx+1;
end

    % if strcmpi(mode_type, 'SH')
    % 
    %     fun = @(sj) det((diag(An)-Ac*B)/10);
    % 
    % elseif strcmpi(mode_type, 'H')
    % 
    %     fun = @(sj) det((diag([A([-n:-1],sj), (2./k_star).*(g_sign*k_star + 1/Bd.*k_star.^3),A([1:n],sj)])-Ac*B)/1);
    % end

end