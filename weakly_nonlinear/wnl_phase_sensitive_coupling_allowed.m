function tf = wnl_phase_sensitive_coupling_allowed(targetMode,sourceMode)
%WNL_PHASE_SENSITIVE_COUPLING_ALLOWED Test bar(A_j)*A_k^2 selection rules.

target = mode_spec(targetMode);
source = mode_spec(sourceMode);
mOut = 2*source.m-target.m;
sOut = wnl_wrap_quasifrequency(2*source.s-target.s);
sError = abs(wnl_wrap_quasifrequency(sOut-target.s));
tf = mOut == target.m && ...
    sError <= 256*eps(max([1,abs(target.s),abs(source.s)]));
end

function spec = mode_spec(value)
if isstruct(value) && isfield(value,'spec')
    spec = value.spec;
else
    spec = value;
end
end
