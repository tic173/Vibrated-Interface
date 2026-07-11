
poolobj = parpool('Processes', 24);
fprintf('Number of workers: %g\n', poolobj.NumWorkers);

clc;
clear; clc; close all;


%%
%%% Constants


 g = 9.81;           % background gravity

 h          =   22/1000;                  % [m]
 R          =   35/1000;

tc         =   sqrt(h/g);



sigma      =   72/1000;    % surface tension coefficient
rho        =   997;  % liquid density
 nu         =   1e-3/rho;        % kinematic viscosity of liquid


%%% Nondimensionalization

 omega      =   11*2*pi;                  % [rad/s]

omega_star =   omega*tc;

C          =  nu/(sqrt(g*h^3));
Bd         =  rho*g*h^2/sigma;


At = 0.9976;
eta = 1.81e-2;

R0 = R/h;



  %%  check roots

  g_sign = -1;

  m = 2;

  n = 10;


  f_count = 0;

  NA = 201;


  Ac_range = linspace(0,0.5,NA);
  f_range = 5:0.1:30;
  f_Ac = nan(length(f_range), length(Ac_range));

  for f = f_range

      % for f = 10.0:0.05:20

      omega      =  f*2*pi;                  % [rad/s]
      omega_star =  omega*tc;

      f_count = f_count+1;

      fun = @(Ac) Ac-Floquet_cylindrical_pinned(Ac,omega_star, R0, m, C, Bd ,At,eta , n ,'SH', g_sign );

       fun1 = @(Ac) Ac-faradayFloquet_RT_boundary_GenEIG_cylindrical_pinned(Ac,omega_star, R0, m, j, C, Bd ,At,eta , n ,'SH', g_sign );

      a_count = 0;

      parfor Ai = 1:NA

          f_Ac(f_count,Ai) = fun(Ac_range(Ai));

      end

      save('Ac_RT_roots','f_Ac','-v7.3');

      % plot3(Ac_range, f*ones(length(Ac_range),1),f_Ac(f_count,:),'bx'); hold on;
      % plot3(Ac_range, f*ones(length(Ac_range),1),zeros(length(Ac_range),1),'k--'); hold on;
      % % ylim([0 5])
      % set(gca,'FontSize',12)
      % drawnow;
      % 
      % % xlabel('freqency (Hz)');
      % % ylabel('Critical acceleration (m^2/s)')
      % 
      % % title(['l = ' num2str(m) ', n = 1'])

  end



%%

  NA = 201;

  Ac_range = linspace(0,0.5,NA);
  f_range = 5:0.1:30;

load('Ac_RT_roots_m3.mat')


  [X,Y]= meshgrid(f_range(1:end),Ac_range);

  figure; colormap(bluered)

  f_plot = f_Ac;
  % f_plot(f_plot < -0.05 | f_plot > 0.05) = NaN;

  pcolor(X,Y,f_plot');
  shading interp;

  caxis([-0.02 0.02])
  %
  xlim([5 30]);
  ylim([0 1])

  xlabel('f');
  ylabel('a')


  Ac_pinned_SH = nan(1,length(f_range));

  for f_count = 1:length(f_range)

      % Ac_pinned_SH(f_count) = f_Ac(f_count,:);

      if any(f_Ac(f_count,1:end-1).*f_Ac(f_count,2:end) < 0)
          disp('Sign change exists.')

       [~,idx]=   min(abs(f_Ac(f_count,:)));
          Ac_pinned_SH(f_count) = Ac_range(idx);
      else
          disp('No sign change.')
      end

  end


           % csvwrite(['Cylinder_Ac_RT_pinned_m' num2str(m) '.csv'],[(5:0.1:30)' Ac_pinned_SH(1:length((5:0.1:30)))']);
