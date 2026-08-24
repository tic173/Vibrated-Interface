function s = wnl_wrap_quasifrequency(s)
%WNL_WRAP_QUASIFREQUENCY Wrap Floquet quasifrequency into (-1/2, 1/2].

s = mod(s + 0.5, 1.0) - 0.5;
tol = 64 * eps(max(1.0, abs(s)));
if s <= -0.5 + tol
    s = 0.5;
end
end
