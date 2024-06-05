function frag = findFragments( fail, sr, hr )
% INPUT VARS: fail - 1 if connected else 0
%             hr - horizons
%             hn - horizon position pointer from 1 to numel(fail)+1
% OUTPUT VARS: frag: the fragments of a solid

frag = cell(0);
deja = [];
for i = 1:1:size(sr, 1)-1
    if(~ismember(i,deja))
        efra0 = [];
        efra1 = i;
        while(numel(efra1)~=numel(efra0))
            %length(efra1)
            emore = efra1(length(efra0)+1:1:length(efra1),1);
            efra0 = efra1;
            for j = 1:1:length(emore)
                k = emore(j);
                maytouched = hr(sr(k):1:sr(k+1)-1,1);
                isconnected = fail(sr(k):sr(k+1)-1,1)==1;
                efra1 = [efra1;maytouched(isconnected,1)];
            end
            efra1 = unique(efra1,'stable');
        end
        deja = [deja;efra1];
        frag = [frag;efra1];
    end
end
end