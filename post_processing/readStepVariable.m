function var = readStepVariable(model, stepNumber, vChar)
% 读入 RepresentativeRadius
fileName = [model.fileDir, model.partName, vChar, num2str(stepNumber),'.bin'];
fid = fopen(fileName, 'rb');
switch lower(vChar)
    case {'displacement', 'damage'}
        [var,~] = fread(fid,'float64');
        var     = reshape(var, model.pn, numel(var)/model.pn);
    case {'fail'}
        [var,~] = fread(fid,'bit8');
        var     = var(:);
end
fclose(fid);
end