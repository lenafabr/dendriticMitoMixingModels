% -------------
% This is an example script for calculating steady-state mitochondrial densities on a
% network, using a model that includes transport, pausing, fusion, fission
% To run: you will need the following dependencies (clone from github):
% https://github.com/lenafabr/dendriticMitoMixing (current package)
% https://github.com/lenafabr/networktools
% Make sure y
% -------------

clear

% ---------- SET THESE TO MATCH YOUR DIRECTORY LOCATIONS --------------
% define the directory where the dendriticMitoMixing package is located
dendriticMitoMixingDir = '~/UCSD/proj/dendriticMitoMixing/dendriticMitoMixingModels_public/';
% define the directory where the networktools package is located
networktoolsDir = '~/UCSD/proj/networktools/';

%% Make matlab aware of  subdirectories and dependencies
cd(dendriticMitoMixingDir)
addpath(genpath('./'))
addpath(networktoolsDir)

%% Load in an HSE neuron dendritic arbor structure, from a .net file
% Treat it as a 2D network
netfile = 'MCFO-HSE-1_UM_exptrad.net';
NT = NetworkObj(['./examples/',netfile],struct('dim',2));

% rearrange edge directions so they always point from upstream node to
% downstream node
% SET rootnode to be the root of your tree
NT.rootnode = 80; 
directedTreeEdges(NT,NT.rootnode);
% this is the proximal trunk of your network
trunkedge = NT.nodeedges(NT.rootnode,1);


% SET whether you want to recalculate radii or use pre-saved ones
% true if you want to use balanced radii (Da Vinci Law, L/D sister area ratios)
% false if you want to use measured radii, with radii.^2 saved in the .net file
set_balanced_radii = true;
% SET the trunk radius in um
rtrunk = sqrt(NT.edgevals(trunkedge,2));

if (set_balanced_radii) % recalculate balanced radii
    % Da Vinci Law: r1^2 + r2^2 = r0^2
    a = 2;
    rm = 0; % minimal allowed radius (in um)

    % calculate length, depth for all subtrees
    [stL,~,stD,~] = setSubtreeInfo(NT,trunkedge,a,'L/D');

    % Split branch areas in proportion to bushiness: r1^2/r2^2 = (L1/D1)/(L2/D2)    
    radii= setRadiiWithRm(NT,trunkedge,2,rm,rtrunk,stL./stD);

else
    % use radius of each edge loaded from the .net file
    radii = sqrt(NT.edgevals(:,2))';
end

% plot the network, with node and index labels accessed by hovering mouse pointer
figure
NT.plotNetwork(struct('datatipindex',1));
title(sprintf('Network structure: %s',netfile),'Interpreter','none')

%% -------------------
% SET UP DYNAMIC PARAMETERS
% describing mitochondrial motion and interaction
% -----------------------

param = struct();

% anterograde mito flux: # units entering trunk per second
param.kp=1/60;


% velocity of anterograde and retrograde moving mitos (um/sec)
param.v=0.5;

% Pausing occurs at rate: As/radius^gammas (units of per sec)
param.gammas = 0.8;
param.As=0.1;

% rate of a paused mitochondrion restarting and becoming mobile again
param.kr = 0.0034;

% Rate of fission (from each side of a cluster); units of per sec
param.kb = 5e-4;
% Probability of fusion when two mitos pass: Au/radius^gamma
param.gamma = 2;
param.Au = 0.1;

% cut off fusion probabilities, so they cannot go above 1
param.probmax = 1;

% reactivation rate (in both directions together)
% for a unit-size cluster left behind after fission
param.ka = 2*param.kb;

% Below are some definitions just for calculating metrics
% distcutoff defines the distal region of the tree (for calculatin distal
% enrichment)
% graph distance to trunk must be more than this fraction of the max
% graph distance
param.distcutoff = 0.7;

% volume of a mitochondrial unit (um^3)
% (for calculating volume fraction, distal enrichment)
param.mitovol = 0.5;

% -------------------------------
%% Compute mitochondrial distribution
% -------------------------------

distinfo = getMitoDensitiesMetrics(NT,radii,param);

Pu0 = min(param.Au/radii(trunkedge)^param.gamma,param.probmax);
kstop0 = param.As/radii(trunkedge)^param.gammas;

disp(sprintf('Proximal fusion probability, stopping rate: %g %g', Pu0,kstop0))
disp(sprintf('average cluster size, volume fraction, distal enrichment: %g %g %g',distinfo.avgclustsize, distinfo.volfrac,distinfo.distenrich))

% -------------------------------
%% Plot volume density of mitochondria on different tree branches
% -------------------------------
% M1vals = # of mito units per length
% voldens = volume of mitos per volume of branch
voldens = distinfo.M1vals*param.mitovol./(pi*radii.^2);

figure
NT.plotNetworkField(radii,voldens);
title('Predicted volume density of mitos')
colorbar