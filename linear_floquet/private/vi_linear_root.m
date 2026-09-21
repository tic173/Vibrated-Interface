function result = vi_linear_root(p,kStar,options,N,extraGuesses)
% Multistart complex root search using the existing Schur coefficients.
% A full complex exponent is searched with harmonics -N:N; no imposed H/SH
% frequency is added to it during reconstruction.
if nargin<5, extraGuesses=[]; end
indices=(-N:N).';
forcing=vi_floquet_acceleration_matrix(numel(indices),p.phase);
% Finite-depth inviscid frequency provides seeds, not the viscous answer.
xi=(1-p.At)/(1+p.At);
frequencySquared=((1-xi)*p.gSign*kStar+kStar^3/p.Bd)/ ...
    (coth(kStar*p.depths(1))+xi*coth(kStar*p.depths(2)));
inviscid=sqrt(complex(-frequencySquared));
damping=2*p.C*(1+p.eta)/(1+xi)*kStar^2;
if isempty(options.initialGuesses)
    span=max(1,p.omegaStar);
    [rr,ii]=ndgrid(span*[-.5 -.05 .05 .5],p.omegaStar*[-.49 -.25 0 .25 .49]);
    guesses=[rr(:)+1i*ii(:); -damping+inviscid; -damping-inviscid];
else
    guesses=options.initialGuesses(:);
end
guesses=unique([extraGuesses(:);guesses]);
settings=optimset('Display','off','MaxIter',options.maxIterations, ...
    'MaxFunEvals',options.maxIterations*8,'TolX',1e-10,'TolFun',1e-10);
candidates=struct('exponent',{},'zeta',{},'residual',{},'edgeEnergy',{},'exitFlag',{});
failures=cell(numel(guesses),1);
% Singular trial vertical systems can occur at isolated internal Stokes
% poles. They are rejected; warnings are restored when this call returns.
w1=warning('off','MATLAB:nearlySingularMatrix');
w2=warning('off','MATLAB:singularMatrix');
cleanup=onCleanup(@() restore_warnings(w1,w2));
for j=1:numel(guesses)
    try
        [xy,~,flag]=fsolve(@residual,[real(guesses(j));imag(guesses(j))],settings);
        s=xy(1)+1i*xy(2);
        if ~isfinite(s), continue; end
        % Fold aliases, then solve again at the canonical exponent so the
        % vector and its finite harmonic window belong to the SAME root.
        canonical=real(s)+1i*wrap(imag(s),p.omegaStar);
        if abs(canonical-s)>1e-7*p.omegaStar
            [xy,~,flag]=fsolve(@residual,[real(canonical);imag(canonical)],settings);
            s=xy(1)+1i*xy(2);
        end
        [M,S]=matrix(s);
        [~,sv,V]=svd(S,'econ'); z=V(:,end);
        r=max(norm(S*z)/max(norm(S,'fro')*norm(z),eps), ...
            norm(M*z)/max(norm(M,'fro')*norm(z),eps));
        r=max(r,sv(end,end)/max(sv(1,1),eps));
        edge=sum(abs(z([1 end])).^2)/sum(abs(z).^2);
        if flag<=0 || ~isfinite(r) || r>options.rootTolerance || ...
                edge>options.edgeTolerance || abs(imag(s))>(.5+1e-6)*p.omegaStar
            failures{j}='Rejected by convergence, matrix residual, alias, or harmonic-edge check.';
            continue
        end
        canonical=real(s)+1i*wrap(imag(s),p.omegaStar);
        duplicate=false;
        for kk=1:numel(candidates)
            old=candidates(kk).exponent;
            if abs(real(old)-real(s))<1e-6*max(1,p.omegaStar) && ...
                    abs(wrap(imag(old)-imag(canonical),p.omegaStar))<1e-6*max(1,p.omegaStar)
                duplicate=true;
                if r<candidates(kk).residual
                    candidates(kk)=record(s,z,r,edge,flag);
                end
                break
            end
        end
        if ~duplicate, candidates(end+1)=record(s,z,r,edge,flag); end %#ok<AGROW>
    catch err
        failures{j}=[err.identifier ': ' err.message];
    end
end
if isempty(candidates)
    error('vi_linear:NoRoot', ...
        ['No acceptable root at kh=%.6g (N=%d). Increase harmonics or broaden ', ...
         'numerics.initialGuesses; inspect parameters and boundary conditions.'],kStar,N);
end
% Deterministic conjugate selection: prefer positive quasi-frequency among
% roots with indistinguishable largest growth rates.
growth=real([candidates.exponent]); best=max(growth);
near=find(abs(growth-best)<1e-7*max(1,p.omegaStar));
[~,pick]=max(imag([candidates(near).exponent])); selected=near(pick);
result=candidates(selected); result.indices=indices; result.harmonics=N;
result.candidates=candidates; result.initialGuesses=guesses;
result.rejectedTrials=failures; result.searchComplete=false;

    function [M,S]=matrix(s)
        a=vi_reduced_cylinder_coefficients(N,s,kStar,p.omegaStar, ...
            p.At,p.eta,p.C,p.Bd,'H',p.gSign,p.depths);
        M=diag(a)-p.Ac*forcing;
        if any(~isfinite(M(:))), error('vi_linear:Trial','Nonfinite trial matrix.'); end
        scale=max(max(abs(M),[],2),1);
        S=M./scale;
    end
    function f=residual(x)
        [~,S]=matrix(x(1)+1i*x(2));
        eigenvalues=eig(S); [~,ix]=min(abs(eigenvalues));
        value=eigenvalues(ix); f=[real(value);imag(value)];
    end
end

function s=record(exponent,zeta,residual,edgeEnergy,exitFlag)
s=struct('exponent',exponent,'zeta',zeta,'residual',residual, ...
    'edgeEnergy',edgeEnergy,'exitFlag',exitFlag);
end
function y=wrap(x,omega)
y=mod(x+omega/2,omega)-omega/2;
end
function restore_warnings(w1,w2)
warning(w1); warning(w2);
end
