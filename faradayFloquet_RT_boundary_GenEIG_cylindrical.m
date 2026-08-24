function [Ac, zetaNeutral, diagnostics] = ...
    faradayFloquet_RT_boundary_GenEIG_cylindrical( ...
    omega_star, R0, m, l, C, Bd, At, eta, n, varargin)
% s_star - growth rate
% omega_star - driving frequence
% beta_star - wave number
% n - order of expansion

mode_type       = 'SH';
g_sign = 1;
forcing_phase = 0;


if nargin >=10
    mode_type = varargin{1};
end

if nargin >=11
    g_sign = varargin{2};
end

if nargin >=12
    forcing_phase = varargin{3};
end

 % 
  beta_star = bessel_derivative_root(m, l)/R0;

  % beta_star = besselzero(m, l+1,1)/R0;

beta_star = beta_star(l);


%%

% mode_type
% 
% g_sign


    if strcmpi(mode_type, 'SH')
        Ln  = 2*(n+1);
        B = vi_floquet_acceleration_matrix(Ln, forcing_phase);

        sj = 0.5*1i*omega_star;
     
    elseif strcmpi(mode_type, 'H')
        Ln  =  2*n+1;
        B = vi_floquet_acceleration_matrix(Ln, forcing_phase);

        sj = 0;

    end

    A0 = vi_reduced_cylinder_coefficients(n, sj, beta_star, ...
        omega_star, At, eta, C, Bd, mode_type, g_sign);


    tol = 1e-6;
    [V, D] = eig(diag(A0), B);
    eigenvalues = diag(D);
    nearlyReal = abs(imag(eigenvalues)) < ...
        tol*max(abs(real(eigenvalues)), 1);
    admissible = find(nearlyReal & isfinite(eigenvalues) & ...
        real(eigenvalues) >= 0);
    if isempty(admissible)
        error('faradayFloquet_RT_boundary:NoPositiveThreshold', ...
            ['No finite nonnegative neutral acceleration was found for ', ...
             'm=%d, l=%d, mode=%s.'], m, l, mode_type);
    end
    [Ac, location] = min(real(eigenvalues(admissible)));
    selectedIndex = admissible(location);
    zetaNeutral = V(:, selectedIndex);
    zetaNeutral = zetaNeutral/max(norm(zetaNeutral), eps);

    neutralMatrix = diag(A0)-Ac*B;
    singularValues = svd(neutralMatrix);
    diagnostics = struct();
    diagnostics.modeType = mode_type;
    diagnostics.forcingPhase = forcing_phase;
    diagnostics.harmonicExponent = sj;
    diagnostics.eigenvalues = eigenvalues;
    diagnostics.selectedEigenvalueIndex = selectedIndex;
    diagnostics.relativeSingularResidual = min(singularValues)/ ...
        max(max(singularValues), eps);
    diagnostics.vectorResidual = norm(neutralMatrix*zetaNeutral)/ ...
        (max(norm(neutralMatrix, 2), eps)*norm(zetaNeutral));


end


% function [An] = FD_RT_coeffs(n,s_star,beta_star,omega_star,At,eta,C,Bd,mode_type,g_sign)
% 
% 
% if strcmpi(mode_type, 'SH')
%     Ln  = -n-1:n;
% 
% elseif strcmpi(mode_type, 'H')
%     Ln  =  -n:n;
% 
% end
% 
% 
% xi = (1-At)/(1+At);
% 
% An = zeros(length(Ln),1);
% 
% n_idx = 1;
% 
% for nj = Ln
% 
%     q1 =  sqrt(1+(s_star+1i*nj*omega_star)/(C*beta_star^2));
%     q2 =  sqrt(1+(s_star+1i*nj*omega_star)*xi/(eta*C*beta_star^2));
% 
% 
% 
%     % reduce to infinite height for large wave number
% 
%     L1 = [1,1,0,0];
%     L2 = [1,1,-1,-1];
%     L3 = [1,q1,1,q2];
%     L4 = [2,(q1^2+1),-2*eta,-eta*(1+q2^2)];
% 
%     LHS = [ L1;L2;L3;L4];
% 
%     RHS = [s_star+1i*nj*omega_star;0;0;0];
% 
%     coeffs0 = (LHS)\RHS;
%     coeffs = [coeffs0(1); 0; coeffs0(2);0; 0; coeffs0(3);0;coeffs0(4)];
% 
% 
%     dwdz =  (coeffs(1)-coeffs(2)+coeffs(3)*q1-coeffs(4)*q1)*beta_star;
% 
%     d3wdz3 =  (coeffs(1)-coeffs(2)+coeffs(3)*q1^3-coeffs(4)*q1^3)*beta_star^3;
% 
% 
%     An(n_idx) = 2*dwdz*((s_star+1i*nj*omega_star)/(beta_star^2)+3*C*(1-eta)*(1+At)/2/At)-2*C*d3wdz3/(beta_star^2)*(1+At)/2/At+2*(g_sign+beta_star^2*(1+At)/Bd/2);
% 
%     if nj ==0 && strcmpi(mode_type, 'H')
%         An(n_idx) = (2./beta_star).*(g_sign*beta_star + (1+At)/2/Bd.*beta_star.^3);
%     end
% 
%     n_idx = n_idx+1;
% end
% 
% 
% end



function [An] = legacy_FD_RT_coeffs_unused(n,s_star,beta_star,omega_star,At,eta,C,Bd,mode_type,g_sign) %#ok<DEFNU>
% Retained only for source-history comparison. The executable solver uses
% vi_reduced_cylinder_coefficients so threshold and growth calculations
% cannot silently diverge.


if strcmpi(mode_type, 'SH')
    Ln  = -n-1:n;

elseif strcmpi(mode_type, 'H')
    Ln  =  -n:n;

end


xi = (1-At)/(1+At);

An = zeros(length(Ln),1);

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


        if  abs(real(exp(beta_star*(-1+q1))))>1e2 && abs(real(exp(beta_star*(-1+q2))))>1e2

            L1 = [1,1,1,0,0,0];
            L2 = [1,1,1,-1,-1,-1];
            L3 = [1,-1,q1,-1,1,q2];
            L4 = [0,0,1,0,0,-eta];

            % L5 = [exp(-beta_star*(2)),1,0,0,0,0];
            % 
            % L6 = [(1+q1)*exp(-beta_star*2),(-1+q1),0,0,0,0];

            L5 = [exp(-beta_star*(1)),-exp(beta_star*(1)),q1*exp(-beta_star*q1),0,0,0];

            L6 = [0,0,0, exp(beta_star*(1)),-exp(beta_star*(-1)),-q2*exp(-beta_star*q2)];

            LHS = [ L1;L2;L3;L4;L5;L6];

            RHS = [s_star+1i*nj*omega_star;0;0;0;0;0];

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

        RHS = [s_star+1i*nj*omega_star;0;0;0];

        coeffs0 = (LHS)\RHS;
        coeffs = [coeffs0(1); 0; coeffs0(2);0; 0; coeffs0(3);0;coeffs0(4)];

    end


    dwdz =  (coeffs(1)-coeffs(2)+coeffs(3)*q1-coeffs(4)*q1)*beta_star;

    d3wdz3 =  (coeffs(1)-coeffs(2)+coeffs(3)*q1^3-coeffs(4)*q1^3)*beta_star^3;

    An(n_idx) = 2*dwdz*((s_star+1i*nj*omega_star)/(beta_star^2)+3*C*(1-eta)*(1+At)/2/At)-2*C*(1-eta)*d3wdz3/(beta_star^2)*(1+At)/2/At+2*(g_sign+beta_star^2*(1+At)/Bd/2/At);

    if nj ==0 && strcmpi(mode_type, 'H')
         An(n_idx) = (2./beta_star).*(g_sign*beta_star + (1+At)/2/Bd.*beta_star.^3);
    end

    n_idx = n_idx+1;
end


end
