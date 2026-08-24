function [s_star,zeta_all,w1,w2,w1p,w2p,diagnostics] = faradayFloquet_RT_GenEIG_cylindrical(Ac, omega_star, R0, m, l, C, Bd, At, eta, n, varargin)
%FARADAYFLOQUET_RT_GENEIG_CYLINDRICAL Linear Floquet exponent and mode.
%
% s_star is the complex Floquet exponent. real(s_star) is the
% dimensionless growth rate and imag(s_star) is the modal angular
% frequency. The optional inputs are
%   mode_type   'H' or 'SH'
%   g_sign      signed gravity parameter
%   s_guess     initial guess for the complex Floquet exponent
%   phase       phase in cos(omega*t+phase), matching the full operator
%
% Velocity fields are reconstructed when 3--6 outputs are requested.
% Calling with seven outputs returns diagnostics and skips the large z-grid
% reconstruction; w1, w2, w1p, and w2p are then empty.

mode_type       = 'SH';
g_sign = 1;
s_guess = [];
forcing_phase = 0;


if nargin >=11
    mode_type = varargin{1};
end

if nargin >=12
    g_sign = varargin{2};
end

if nargin >=13
    s_guess = varargin{3};
end

if nargin >=14
    forcing_phase = varargin{4};
end


beta_star = bessel_derivative_root(m, l)/R0;

beta_star = beta_star(l);


%%  Construct the matrix A


[A0] = @(sj) vi_reduced_cylinder_coefficients(n, sj, beta_star, ...
    omega_star, At, eta, C, Bd, mode_type, g_sign);

if strcmpi(mode_type, 'SH')
    Ln  = 2*(n+1);
    B = vi_floquet_acceleration_matrix(Ln, forcing_phase);

elseif strcmpi(mode_type, 'H')
    Ln  =  2*n+1;
    B = vi_floquet_acceleration_matrix(Ln, forcing_phase);

end


%%    Solve

fun = @(sj) det((diag(A0(sj))-Ac*B)/20);

opt=optimset('Maxiter',2000,'TolX',1e-10,'Tolfun',1e-10);


  % s_star = fsolve(@(sj) fun(sj),  3+1i*1.6365 );

% s_star = fsolve(@(sj) fun(sj),  1+1*beta_star/2.5-1*beta_star^2*C/2+sqrt(beta_star+C*(beta_star).^3)*1i, opt);

 
if isempty(s_guess)
    s_guess = 1+beta_star/3-C + ...
        0.4i*sqrt(beta_star+C*beta_star^3);
end
validateattributes(s_guess, {'numeric'}, {'scalar', 'finite'});
[s_star, determinantResidual, exitFlag, solveOutput] = ...
    fsolve(@(sj) fun(sj), s_guess, opt);

  % s_star = fsolve(@(sj) fun(sj), 1+1*beta_star/3-1*beta_star^2*C*1+sqrt(beta_star+C*(beta_star).^3*1)*1i*0.1, opt);

 % s_star = fsolve(@(sj) fun(sj), 1, opt);   % m=0


% if imag(s_star)<1e-4 && real(s_star)<0
%     s_star = fsolve(@(sj) fun(sj), abs(s_star), opt);
% end


%%   Periodic Harmonic Components

modeMatrix = diag(A0(s_star))-Ac*B;
singularValues = svd(modeMatrix);
relativeSingularResidual = min(singularValues) / ...
    max(max(singularValues), eps);
zeta_all = null(modeMatrix);
if isempty(zeta_all)
    [~, ~, rightVectors] = svd(modeMatrix, 'econ');
    zeta_all = rightVectors(:, end);
elseif size(zeta_all, 2) > 1
    zeta_all = zeta_all(:, 1);
end

w1 = [];
w2 = [];
w1p = [];
w2p = [];
diagnostics = struct();
diagnostics.initialGuess = s_guess;
diagnostics.determinantResidual = abs(determinantResidual);
diagnostics.exitFlag = exitFlag;
diagnostics.solveOutput = solveOutput;
diagnostics.relativeSingularResidual = relativeSingularResidual;
diagnostics.singularValues = singularValues;
diagnostics.acceleration = Ac;
diagnostics.forcingPhase = forcing_phase;

if nargout >= 3 && nargout <= 6
    [w1,w2,w1p,w2p] = FD_RT_BCcoeffs(n,s_star,beta_star, ...
        omega_star,At,eta,C,Bd,mode_type,g_sign,zeta_all);
end

 % w1 = 0; w2 =0; w1p=0;w2p=0; zeta_all=0;

% plot(t/(2*pi/omega_star),real(exp(growthrate_FD(j)*t).*transpose(sum(zeta_all(j,2:end).*exp(1i*omega_star*(-n:n).*t'),2)/sum(zeta_all(j,2:end).*exp(1i*omega_star*(-n:n).*(0+1*phase)),2))),'linewidth',2,'LineStyle','-','Color','r');


end


function [An] = legacy_FD_RT_coeffs_unused(n,s_star,beta_star,omega_star,At,eta,C,Bd,mode_type,g_sign) %#ok<DEFNU>
% Retained only for source-history comparison. The executable solver uses
% vi_reduced_cylinder_coefficients.


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

    % % reduce to infinite height for large wave number
    %
    % L1 = [1,1,0,0];
    % L2 = [1,1,-1,-1];
    % L3 = [1,q1,1,q2];
    % L4 = [2,(q1^2+1),-2*eta,-eta*(1+q2^2)];
    % 
    % LHS = [ L1;L2;L3;L4];
    % 
    % RHS = [s_star+1i*nj*omega_star;0;0;0];
    % 
    % coeffs0 = (LHS)\RHS;
    % coeffs = [coeffs0(1); 0; coeffs0(2);0; 0; coeffs0(3);0;coeffs0(4)];


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


     dwdz =  (coeffs(1)-coeffs(2)+coeffs(3)*q1-coeffs(4)*q1)*beta_star;

    % dwdz =  (coeffs(5)-coeffs(6)+coeffs(7)*q2-coeffs(8)*q2)*beta_star;

    d3wdz3_1 =  (coeffs(1)-coeffs(2)+coeffs(3)*q1^3-coeffs(4)*q1^3)*beta_star^1;

    d3wdz3_2 =  (coeffs(5)-coeffs(6)+coeffs(7)*q2^3-coeffs(8)*q2^3)*beta_star^1;

     % An(n_idx) = 2*dwdz*((s_star+1i*nj*omega_star)/(beta_star^2)+3*C*(1-eta)*(1+At)/2/At)-2*C*d3wdz3/(beta_star^2)*(1+At)/2/At+2*(g_sign+beta_star^2*(1+At)/Bd/2);

       % An(n_idx) = 2*dwdz*((s_star+1i*nj*omega_star)/(1*beta_star^2)+3*C*(1-eta)*(1+At)/2/At)-2*C*(1-eta)*d3wdz3_1/(beta_star^0)*(1+At)/2/At+2*(g_sign+beta_star^2*(1+At)/Bd/2/At);

     An(n_idx) =      2*dwdz*( (s_star+1i*nj*omega_star)/beta_star^2 ...
    + 3*C*(1-eta)*(1+At)/(2*At) ) ...
    - C*(1+At)/At*( d3wdz3_1 - eta*d3wdz3_2 ) ...
    + 2*(g_sign + beta_star^2*(1+At)/(2*Bd*At));

      % An(n_idx) =  An(n_idx) *10;

    if nj ==0 && strcmpi(mode_type, 'H')
        An(n_idx) = (2./beta_star).*(g_sign*beta_star + (1+At)/2/At/Bd.*beta_star.^3);
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
