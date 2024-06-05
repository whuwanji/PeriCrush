function [frag, fragVol] = volumeDistribution(model, stepNumber, rateDx)
% 根据时间步确定碎片所包含物质点编号以及碎片的体积
pv = model.ParticleVolume;
hr = model.Horizon;
% 加载所有近场域的质点数目
hn = model.HorizonParticleNumber;
fail      = readStepVariable(model, stepNumber, 'fail');        % 读入第stepNumber步的损伤
dis      = readStepVariable(model, stepNumber, 'displacement');        % 读入第stepNumber步的损伤
dx = model.dx;
sr = zeros(model.pn+1,1);
sr(1) = 1;
for i = 1:1:model.pn
    sr(i+1) = sr(i) + hn(i);
end
hi = zeros(size(hr));
for i = 1:1:numel(sr)-1
    hi(sr(i):1:sr(i+1)-1,1) = i;
end
coor = model.Coordinate + dis;
br = sqrt(sum((coor(hr,:)-coor(hi,:)).^2, 2))<(rateDx*dx);
fail = fail==1&br;
frag = findFragments( fail, sr, hr );
loadp = [];
issingle = false(size(frag));
for i = 1:1:numel(frag)
    if(length(frag{i})<5)
        issingle(i)=true;
    end
end
frag(issingle)=[];
for i = 1:1:numel(frag)
    efra = setdiff(frag{i},loadp);
    frag{i} = efra;
end
fragVol = zeros(size(frag));
for i = 1:1:numel(frag)
    efra = frag{i};
    fragVol(i) = sum(pv(efra));
end
[fragVol, idx]=sort(fragVol,'descend');
frag = frag(idx);
end