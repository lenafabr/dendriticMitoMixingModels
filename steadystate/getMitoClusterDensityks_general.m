function [M0vals, M1vals,rhoSdvals,rhoSvals] = getMitoClusterDensityks_general(rhoWvals,Pu,alpj,deltaj,ksj,kr,v,Nmax)

% get linear density of mito mass and linear density of clusters on each edge
% includes arrest of motile mitochondria  
% after having set motile mito density with setMotileMitoConcFromRadii.m
% and after having set the alpj = v*Au/(2*kb*r^gamma) values
% rhoSvals has the densities of stationary clusters of different sizes up
% to Nmax
% ksj array of stop rate per branch ksj=As./radii(j)^gammas
% kb rate of breaking
  
% convert everything to row vectors
if (size(alpj,1)>1)
    alpj = alpj';
end
if (size(alpj,1)>1)
    deltaj = deltaj';
end
if (size(rhoWvals,1)>1)
    rhoWvals = rhoWvals';
end
if (size(ksj,1)>1)
    ksj = ksj';
end
if (~exist('Nmax','var'))
    Nmax = 20; % default max cluster size
end
rhoSdvals=ksj.*rhoWvals./(kr+v*Pu.*rhoWvals);
% linear density of mito mass on each edge
alrho = alpj.*rhoWvals;
deltarho=deltaj.*rhoWvals;
badind = find(alrho>1);
M1vals = rhoWvals/2 + ((rhoWvals/2+rhoSdvals).*(1 + deltarho))./((alrho-1).^2);
M1vals(badind) = inf;

% linear density of all clusters on each edge
M0vals = rhoWvals/2 + ((rhoWvals/2+rhoSdvals).*(1+ deltarho))./(1-alrho);
M0vals(badind) = inf;

if (nargout>2)
    % clusters of different sizes
    nvals = (1:Nmax)';
    arho = alpj.*rhoWvals;
    drho=deltaj.*rhoWvals;
    rhoSvals = (rhoWvals./2+rhoSdvals).*(1+drho).*(arho).^(nvals-1);
    rhoSvals(1,:) = deltaj.*rhoWvals.^2/2 + deltaj.*rhoWvals.*rhoSdvals;
end

end

