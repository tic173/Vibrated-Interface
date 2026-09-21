function result = vi_linear_floquet(cfg)
%VI_LINEAR_FLOQUET Compute growth and dynamics for every requested mode.
% For a scan followed by ONLY the fastest mode, use vi_linear_most_unstable.
result=vi_linear_dynamics(vi_linear_growth(cfg));
end
