function  [output] = getMitoDensitiesMetrics(NT,radii,param)
% This function computes steady-state mitochondrial densities on all
% network branches
% And also some metrics, such as the overall volume fraction for
% mitochondria, the average cluster size, and the distal enrichment
% [Historic note: copied from getAllMetricsFuseArrest_gen_param.m]


% Input: 
% NT = network object structure (defined in https://github.com/lenafabr/networktools)
% radii = 1 x NT.nedge list, radius for each edge of network
% param = structure defining various parameters for model

% Output: structure containing the following fields
% rhoWvals = 1 x nedge array; linear density of motile mitos, on each edge
% rhoSdvals = 1 x nedge array; linear density of paused mitos, on each edge
% rhoSvals = Nmax x nedge; linear density of stopped memoryless clusters of each size, on each edge
% avgclustsize = average cluster size
% volfrac = mito volume fraction
% distenrich = ratio of volume fraction in distal region vs proximal trunk
% M0vals = total linear density of mito clusters on each edge
% M1vals = total linear density of mito units on each edge
% distedges = which edges are considered distal

%% legacy check for compatibility with old code
if (isfield(param,'kfiss'))
    error('Please update your code to use ka instead of kfiss for the rate of reactivation of a unit left after fission')
end

%% Set up default input parameters
% the defaults are overwritten by values from param whenever supplied by
% user
gamma = 2; % radius scaling exponent for fusion
gammas = 0.8; % radius scaling exponent for pausing
Au = 0.9; % fusion rate prefactor
As = 0.09; % pausing rate prefactor
kb = 0.003; % fission rate for clusters of size>1
kp = 1/60; % production rate (# mitos entering per sec)
ka = 2*kb; % restarting rate of size one mitos left behind by fission
kr = 0.001; % restarting rate for paused mitos
v = 0.4; % mito velocities
mitovol = 0.5; % mito unit volume

lencorrect = false; % do finite length correction for trunk mobile density?
Aubuffer = false; % no fusion at distalmost and proximalmost edges

distcutoff = 0.7; %cutoff defining distally enriched branches (as fraction of max branch distance from root)
Nmax = 20; % max cluster size to consider
% cut off probabilities Pu and Puk at a maximum value
% by default do not cut off
probmax = inf;

% stopping rate on each edge (can be set explicitly as a vector in param, or use As,
% gammas)
clear ksj

% unpack the input param structure, overwriting defaults as needed
v2struct(param)

if (size(radii,1)>size(radii,2))
    % convert radii to row vector if it was supplied as a column
    radii = radii';
end

%%
% compute densities of walking (mobile) mitochondria on each edge
% rhoWvals describes total density of anterograde + retrograde mitos
trunkedge = NT.nodeedges(NT.rootnode,1);
if (lencorrect)
    % do finite length correction
    rhoWtrunk = getMotileDensityTrunk_lencorr(param,radii(trunkedge));
else
    rhoWtrunk = 2*kp/v;
end
[rhoWvals] = setMotileMitoConcFromRadii(NT,trunkedge,rhoWtrunk,radii);

%% probability of fusing into larger clusters (on each network edge)
Pu=(Au./radii.^gamma);
% optionally, cut off probabilities at a certain maximum value
Pu=min(Pu,probmax);


% alpha parameter on each edge
alpj = v*Pu./(2*kb);
% delta parameter on each edge
deltaj=v*Pu./(ka);
% stopping rate on each edge (if not already set)
if (~exist('ksj','var'))    
    ksj = As./(radii.^gammas);
end

% get linear density of clusters and total mass on each edge
% rhoSdvals = densities of mitos in dynamic stationary states (paused,
% remembering direction of motion)
% rhoSvals = densities of mitos in fully stationary clusters of different
% size
% M0vals = total linear density of clusters
% M1vals = total linear density of mitochondrial unit count (mass)
if (~exist('Nmax','var'))
    Nmax = 20;
end
[M0vals, M1vals,rhoSdvals,rhoSvals] = getMitoClusterDensityks_general(rhoWvals,Pu,alpj,deltaj,ksj,kr,v,Nmax);

%% average cluster size
if (isinf(max(M1vals)))
  avgclustsize = 1e10;
else
    avgclustsize = sum(M1vals.*(NT.edgelens'))/sum(M0vals.*(NT.edgelens'));
end
%% overall mito volume fraction
volfrac = (mitovol*M1vals*NT.edgelens)/(pi*radii.^2*NT.edgelens);

%% distal enrichment
lend = pathLengthEdges(NT);
distedges = find(lend>distcutoff*max(lend));
% distal avg volume density
Vdist = sum(NT.edgelens(distedges).*radii(distedges).^2);
cdist = sum(M1vals(distedges).*NT.edgelens(distedges))/Vdist;
% primary trunk avg volume density
cprox = M1vals(trunkedge)/(radii(trunkedge).^2);

distenrich = cdist/cprox;


%% pack outputs into structure
output = v2struct(avgclustsize,volfrac,distenrich,M0vals,M1vals,rhoWvals,rhoSdvals,rhoSvals,distedges);

end
