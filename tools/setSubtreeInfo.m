function [stL,stEta,stD,muvals,stN] = setSubtreeInfo(NT,trunkedge,a,splittype,stL,stEta,stD,muvals,stN)
% walk down a tree, starting with the trunkedge
% compute information on subtrees based on tree topologies 
% no initial assumptions made on the absolute radii. No assumptions on
% mitochondrial transport.
% does assume a scaling exponent alpha
% stL = total length of subtree below each junction
% stEta = eta value for the subtree below each junction
% defined by: V_subtree = eta*r_trunk^2
% muvals = r_1^alpha/r_0^alpha (daughter vs parent) for each junction
% stN = number of distal tips in each subtree

% if input stL, stEta, muvals arrays not provided (or empty), initialize to zeros
if (~exist('stL','var'))
    % initialize all arrays to 0
    stL = zeros(1,NT.nedge);    
end
if (~exist('stEta','var'))
    % initialize all arrays to 0
    stEta = zeros(1,NT.nedge);    
end
if (~exist('stD','var'))
    % initialize all arrays to 0
    stD = zeros(1,NT.nedge);    
end
if (~exist('muvals','var'))
    % initialize all arrays to 0
    muvals = zeros(1,NT.nedge);    
end
if (~exist('stN','var'))
    % initialize all arrays to 0
    stN = zeros(1,NT.nedge);    
end
if (isempty(stL))
    % initialize all arrays to 0
    stL = zeros(1,NT.nedge);    
end
if (isempty(stEta))
    % initialize all arrays to 0
    stEta = zeros(1,NT.nedge);    
end
if (isempty(stD))
    % initialize all arrays to 0
    stD = zeros(1,NT.nedge);    
end
if (isempty(muvals))
    % initialize all arrays to 0
    muvals = zeros(1,NT.nedge);    
end

%trunk length
ell0 = NT.edgelens(trunkedge);

%% junction node below this edge
junc = NT.edgenodes(trunkedge,2);

if (NT.degrees(junc)==1)
    % reached a terminal node
    % this edge is a single-trunk subtree
    stEta(trunkedge) = ell0;   
    stD(trunkedge) = ell0;
    muvals(trunkedge) = 0;
    stL(trunkedge) = ell0;
    stN(trunkedge) = 1;

elseif (NT.degrees(junc)==2)
    error('not set up to deal with deg 2 nodes yet')
    
elseif (NT.degrees(junc)==3)
   
    edge1 = NT.nodeedges(junc,2);
    edge2 = NT.nodeedges(junc,3);
    
    % compute subtree properties for the daughter trunks
    % including radii splitting (mu) for all downstream junctions
    [stL,stEta,stD,muvals,stN] = setSubtreeInfo(NT,edge1,a,splittype,stL,stEta,stD,muvals,stN);
    [stL,stEta,stD,muvals,stN] = setSubtreeInfo(NT,edge2,a,splittype,stL,stEta,stD,muvals,stN);
    
    % tree length
    stL(trunkedge) = stL(edge1)+stL(edge2)+ell0;    
    
    % number of tips
    stN(trunkedge) = stN(edge1)+stN(edge2);

    % solve for the radii splitting (mu at the current junction)                
    switch(splittype)
        case('LV')
            % split such that L_1/L2 = V_1/V_2 (subtree length proportional
            % to volume)
            muvals(trunkedge) = 1/(1+((stL(edge2)/stEta(edge2))/(stL(edge1)/stEta(edge1)))^(a/2));                       
        case('L/D')
            % split by bushiness, in proportion to L/D
            muvals(trunkedge) = 1/(1+((stL(edge2)/stD(edge2))/(stL(edge1)/stD(edge1)))^(a/2));
        case('L')
            % split by bushiness, in proportion to L/D
            muvals(trunkedge) = 1/(1+((stL(edge2))/(stL(edge1)))^(a/2)); 
        case('equal')
            % split to equal radii
            muvals(trunkedge) = 0.5;
        otherwise
            error(sprintf('Not a valid value for splittype: %s. Must be LV,L/D,L,equal',splittype))
    end
    
    % tree volume (relative to trunk area)
    stEta(trunkedge) = ell0 + stEta(edge1)*muvals(trunkedge)^(2/a) + stEta(edge2)*(1-muvals(trunkedge))^(2/a);
    
    % subtree depth (same as eta for alpha=2);
    stD(trunkedge) = ell0 + (stL(edge1) + stL(edge2))/(stL(edge1)/stD(edge1) + stL(edge2)/stD(edge2));
else
    error('not set up to deal with node degrees > 3')
end



end