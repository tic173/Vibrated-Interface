function specOut = wnl_combine_spec(model, specs, signs, label)
%WNL_COMBINE_SPEC Return the Floquet block occupied by a nonlinear product.

if ~iscell(specs)
    specs = num2cell(specs);
end
signs = signs(:);
if numel(specs) ~= numel(signs)
    error('wnl_combine_spec:BadSigns', ...
        'There must be one sign for every input specification.');
end

mOut = 0;
sOut = 0;
lambdaOut = 0;
for j = 1:numel(specs)
    mOut = mOut + signs(j) * specs{j}.m;
    sOut = sOut + signs(j) * specs{j}.s;
    lambdaOut = lambdaOut + signs(j) * ...
        wnl_spec_lambda(specs{j});
end
sOut = wnl_wrap_quasifrequency(sOut);
specOut = model.makeSpec(mOut, sOut, label);
specOut.lambda = lambdaOut;
end
