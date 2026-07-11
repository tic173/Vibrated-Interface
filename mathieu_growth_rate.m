function lambda = mathieu_growth_rate(omega0, G, Ac, Omega_star, zeta)

    % One forcing period in nondimensional time
    T = 2*pi/Omega_star;

    % First-order Mathieu system
    rhs = @(t,y) [
        y(2);
        -(omega0^2 + Ac*G*cos(Omega_star*t))*y(1) - 2*zeta*y(2)
    ];

    opts = odeset('RelTol',1e-10,'AbsTol',1e-12);

    % Integrate first basis vector
    [~,Y1] = ode45(rhs,[0 T],[1;0],opts);

    % Integrate second basis vector
    [~,Y2] = ode45(rhs,[0 T],[0;1],opts);

    % Monodromy matrix
    Phi = [Y1(end,:).' Y2(end,:).'];

    % Floquet multipliers
    mu = eig(Phi);

    % Floquet growth rate
    lambda = log(max(abs(mu)))/T;

end