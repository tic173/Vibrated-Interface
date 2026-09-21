function [grid,basis] = vi_linear_spatial(cfg,modes)
g=cfg.geometry; s=cfg.sampling; basis=cell(1,numel(modes));
switch g.type
    case 'cartesian2d'
        x=linspace(0,g.Lx,s.Nx); grid=struct('X',x,'Y',zeros(size(x)));
        for i=1:numel(modes), basis{i}=exp(1i*modes(i).kx*x); end
    case 'cartesian3d'
        [X,Y]=meshgrid(linspace(0,g.Lx,s.Nx),linspace(0,g.Ly,s.Ny));
        grid=struct('X',X,'Y',Y);
        for i=1:numel(modes), basis{i}=exp(1i*(modes(i).kx*X+modes(i).ky*Y)); end
    case 'cylindrical3d'
        [R,theta]=meshgrid(linspace(0,g.radius,s.Nr),linspace(0,2*pi,s.Ntheta));
        grid=struct('X',R.*cos(theta),'Y',R.*sin(theta),'r',R,'theta',theta);
        for i=1:numel(modes)
            m=modes(i).m; k=modes(i).k;
            % Find the continuous radial maximum at endpoints and derivative
            % zeros; the normalization does not depend on the output grid.
            roots=bessel_derivative_root(m,modes(i).radialIndex);
            extremal=besselj(m,[0 roots]);
            normalization=max(abs(extremal));
            basis{i}=besselj(m,k*R).*exp(1i*m*theta)/normalization;
        end
end
end
