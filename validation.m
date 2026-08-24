
%%
%%% Constants

 g = 9.81;

 h          =   22/1000;                  % [m]
 R          =   35/1000;

tc         =   sqrt(h/g);

sigma      =   72/1000;    % surface tension coefficient
rho        =   997;  % liquid density

  nu         =   1e-3/rho;        % kinematic viscosity of liquid

 % nu         =   1e-8/rho;        % kinematic viscosity of liquid

%%% Nondimensionalization
              % [rad/s]

 omega      =   60*2*pi;                  % [rad/s]
omega_star =   omega*tc;

C          =  nu/(sqrt(g*h^3));
Bd         =  rho*g*h^2/sigma;


   At = 0.9976;
   eta = 1.81e-2;

R0 = R/h;

   % At =.5; 
   % eta = 0.3;



%% Theoretical dispersion relation


alpha1 = (1 - At)/2;
alpha2 = (1 + At)/2;

Gamma = alpha2/Bd;

C2 = C;
C1 = eta*(alpha2/alpha1)*C2;   % eta = mu1/mu2

M2 = alpha2*C2;
M1 = eta*M2;

DeltaC = alpha1*C1 - alpha2*C2;

root_theory       = zeros(10,4);
growthrate_RT_theory = zeros(10,4);
growthrate_Faraday_theory = zeros(10,4);

figure;

Ac = 4;

for m = 0:3

    root = bessel_derivative_root(m, 100);


    opt=optimset('Maxiter',2000,'TolX',1e-10,'Tolfun',1e-10);

    for j = 1:10

        k = root(j)*1/R0;

     DeltaC = alpha1*C1 - alpha2*C2;

D_RT = At*k - Gamma*k^3;

Q1 = @(S) sqrt(k^2 + S/C1);
Q2 = @(S) sqrt(k^2 + S/C2);


F113 = @(S) ...
    -(1 + k./S.^2 .* (alpha1 - alpha2 + Gamma*k^2)) .* ...
    (alpha2*Q1(S) + alpha1*Q2(S) - k) ...
    - 4*k*alpha1*alpha2 ...
    + (4*k^2./S).*DeltaC .* ...
    ((alpha2*Q1(S) - alpha1*Q2(S)) + k*(alpha1-alpha2)) ...
    + (4*k^3./S.^2).*DeltaC.^2 .* ...
    (Q1(S) - k).*(Q2(S) - k);

        growthrate_RT_theory(j,m+1) = fsolve(F113,sqrt(D_RT),opt);

         Dk = alpha2*coth(k) + alpha1*coth(k);
                % Dk = alpha2*coth(k) + alpha1;

        omega_n_star = sqrt((At*k + Gamma*k^3)/Dk);
        Gmn = At*k/Dk;
        zeta = 2*(M1 + M2)*k^2;
        growthrate_Faraday_theory(j,m+1) = ...
            mathieu_growth_rate(omega_n_star,Gmn,Ac,omega_star,zeta);

    end
    root_theory(:,m+1) = root(1:10)*h/R;

      % plot(root(1:10)*h/R,real(growthrate_RT_theory(:,m+1)),'o-'); hold on;
      plot(root(1:10)*h/R,real(growthrate_Faraday_theory(:,m+1)),'o-'); hold on;
    drawnow;

end

 % 
   % csvwrite('Faraday_dispersion_case2.csv',[root_theory real(growthrate_Faraday_theory)]);
   % 
   % csvwrite('RT_dispersion_case2.csv',[root_theory real(growthrate_RT_theory)]);
   % 

%%  Faraday growthrate


Nj = 10;

Nm = 4;

Ac=4;
g_sign = -1;


n =10;

growthrate_FD = nan(Nm,Nj);
A_FD = nan(Nm,Nj);

zeta_all = zeros(Nm,Nj,2*n+2);

coeffs_all = zeros(Nm,Nj,2*n+2);


% figure;

for m = 0:Nm-1

      root = bessel_derivative_root(m, 100);

    for j = 1:Nj

      
        % k_star = k_range_FD(j);

           % [growthrate_FD(m+1,j),zeta_all(m+1,j,:), coeffs_all(m+1,j,:)] = faradayFloquet_RT_GenEIG_cylindrical(Ac, omega_star, R0, m, j, C, Bd ,At,eta , n ,'SH', g_sign );

           [growthrate_FD(j,m+1),~, ~] = faradayFloquet_RT_GenEIG_cylindrical(Ac, omega_star, R0, m, j, C, Bd ,At,eta , n ,'SH', g_sign );

    end

    plot(root(1:10)*h/R, real(growthrate_FD(:,m+1)),'linewidth',2,'LineStyle','-'); hold on; drawnow;
    ylabel("Growth rate")

    xlabel("Bessel mode index");

    xlim([0 20])
end


fontsize(gca,16,"points");

title(['acceleration = ' num2str(Ac) 'g, frequency = ' num2str(omega/2/pi) 'Hz' ])

ylim([-inf inf])


legend('$m=0$','$m=1$','$m=2$','$m=3$','$m=4$','$m=5$','Interpreter','latex')


 % csvwrite('RT_dispersion_Floquet_case2.csv',[root_theory real(growthrate_FD)]);

