function [Ac] = faradayFloquet_GenEIG(s_star, omega_star, k_star, h, n, varagin)
% s - Specified growth rate
% omega - driving frequence
% k - wave number (2pi/lambda)
% h - interface height (0 for infinite depth)
% n - order of expansion

%%% Constants
g = 9.81;           % background gravity

% nu = 2e-5;        % kinematic viscosity of liquid

 nu = 2.0e-5;        % kinematic viscosity of liquid
sigma = 2.06e-2;    % surface tension coefficient
rho = 950;          % liquid density

% nu = 1.02e-4;        % kinematic viscosity of liquid
% sigma = 6.76e-2;    % surface tension coefficient
% rho = 950;          % liquid density


%%% Nondimensionalization

tc         =  sqrt(h/g);
C          =  nu/(sqrt(g*h^3));
Bd         =  rho*g*h^2/sigma;
% k_star     =  k.*h;
% s_star     =  s*tc;
% omega_star =  omega*tc;

mode_type       = 'SH';

if nargin ==6
    mode_type = varagin;
end

% if strcmpi(mode_type, 'SH')
%     alpha = 1/2;
% elseif strcmpi(mode_type, 'H')
%     alpha = 0;
% end

 alpha = 0;

tol = 1e-6;

%%

Ls  = length(s_star);
Lk  = length(k_star);


for j = 1:Ls

    sj = s_star(j);

    %%% Linear Stability

    qnh = @(n) sqrt(k_star.^2 + (sj + 1i .* (alpha + n) .* omega_star)/C );


    %%% Infinite depth
    if (h == 0)
        A = @(n) (2./k_star).*(k_star + 1/Bd.*k_star.^3 + C^2.*(qnh(n).^4 + 2.*qnh(n).^2.*k_star.^2 - 4.*qnh(n).*k_star.^3 + k_star.^4));
    else
        %%% Finite depth
        Cn  =  @(n) qnh(n).*(qnh(n).^4 + 2.*qnh(n).^2.*k_star^2 + 5.*k_star.^4);
        Dn  =  @(n) k_star.*(qnh(n).^4 + 6.*qnh(n).^2.*k_star.^2 + k_star.^4);
        num =  @(n) 4.*qnh(n).*k_star.^2.*(qnh(n).^2 + k_star.^2) - Cn(n).*cosh(qnh(n)).*cosh(k_star) + Dn(n).*sinh(qnh(n)).*sinh(k_star);
        dem =  @(n) qnh(n).*cosh(qnh(n)).*sinh(k_star) - k_star.*sinh(qnh(n)).*cosh(k_star);
        A   =  @(n) (2./k_star).*(k_star + 1/Bd.*k_star^3 - C.^2.*(num(n)./(dem(n))));
    end

    if strcmpi(mode_type, 'SH')
        Ln  = 2*(n+1);
        B   =  full( spdiags([ones(Ln,1) zeros(Ln,1) ones(Ln,1)],-1:1,Ln,Ln));
        Aj  =  diag(A([-n-1:n]));
    elseif strcmpi(mode_type, 'H')
        Ln  =  2*n+1;
        B   =  full( spdiags([ones(Ln,1) zeros(Ln,1) ones(Ln,1)],-1:1,Ln,Ln));
        Aj  =  diag(A([-n:n]));
        Aj(n+1,n+1) = (2./k_star).*(k_star + 1/Bd.*k_star.^3);
    end

    [V,D]   =  eig(Aj,B);
    D1       =  (diag(D));

     D       =  real(D1(abs(imag(D1))<tol*abs(real(D1))));

    D       =  sort(D1,'ascend','ComparisonMethod','real')

     Ac(:,j) =  min(D(real(D)>0));

    % Ac(:,j) =  min(D);


    % D       =  D1(real(D1)>0);
    % D       =  sort(1./D,'descend','ComparisonMethod','real');
    % Ac(:,j) =  real(1./D);

end


end