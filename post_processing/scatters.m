function scatters(r, v, markersize, numOfjet)
if(nargin<3)
    markersize = 15;
end
if(nargin<4)
    numOfjet = 11;
end
ndim = size(r,2);
if(ndim==2)
    scatter(r(:,1), r(:,2), markersize, v, 'filled')
    axis equal
    colormap(jet(numOfjet))
    colorbar
elseif(ndim==3)
    scatter3(r(:,1), r(:,2), r(:,3), markersize, v, 'filled')
    axis equal
    colormap(jet(numOfjet))
    colorbar
end
end