function filename = vi_wnl_output_filename(input)
%VI_WNL_OUTPUT_FILENAME Describe a full-cylinder WNL case in its filename.
%
% The filename is deterministic for the operating point and retained modes
% so recovery, coefficient, and initial-condition postprocessing runs resolve
% to the same reusable cache. The complete user and effective inputs,
% including initial amplitudes and phases, remain authoritative inside the
% MAT file.

required = {'dimensional','forcing','numberOfModes','modes'};
for fieldIndex = 1:numel(required)
    if ~isfield(input,required{fieldIndex})
        error('vi_wnl_output_filename:MissingInput', ...
            'input.%s is required.',required{fieldIndex});
    end
end

numberOfModes = input.numberOfModes;
validateattributes(numberOfModes,{'numeric'}, ...
    {'scalar','integer','positive','finite'});
if numel(input.modes) < numberOfModes
    error('vi_wnl_output_filename:ModeCount', ...
        'input.modes contains fewer than input.numberOfModes entries.');
end

modeTokens = cell(1,numberOfModes);
for modeIndex = 1:numberOfModes
    mode = input.modes(modeIndex);
    modeTokens{modeIndex} = sprintf('m%sl%s', ...
        number_token(mode.m),number_token(mode.radialIndex));
end

parts = { ...
    'vi_wnl', ...
    ['ag0-',number_token(input.forcing.analysisAmplitude)], ...
    ['fHz-',number_token(input.dimensional.frequencyHz)], ...
    ['modes-',strjoin(modeTokens,'-')]};
filename = [strjoin(parts,'_'),'.mat'];

if strlength(string(filename)) > 240
    error('vi_wnl_output_filename:FilenameTooLong', ...
        'Automatic output filename exceeds 240 characters: %s',filename);
end
end

function token = number_token(value)
validateattributes(value,{'numeric'}, ...
    {'scalar','real','finite'});
token = lower(sprintf('%.8g',value));
token = strrep(token,'+','');
token = strrep(token,'-','m');
token = strrep(token,'.','p');
end
