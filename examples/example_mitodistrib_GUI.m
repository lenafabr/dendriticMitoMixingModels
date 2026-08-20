% -------------
% This is an example script that calls a GUI which lets you twiddle parameters
% and  calculate steady-state mitochondrial densities on a
% network, using a model that includes transport, pausing, fusion, fission
% To run: you will need the following dependencies (clone from github):
% https://github.com/lenafabr/dendriticMitoMixingModels (current package)
% https://github.com/lenafabr/networktools
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

%% run the GUI
app = mitoDensitiesGUI_exported('NT',NT)

