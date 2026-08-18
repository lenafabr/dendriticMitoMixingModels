function [radii,stV] = setRadiiWithRm(NT,trunkedge,a,rm,rtrunk,splitvals,radii,stV)
% walk down a tree, starting with the trunkedge
% given a radius value for the trunk edge (rtrunk),
% set the radii of the rest of the branches.
% assumes r^2 ratios are proportional to splitvals (eg: splitvals might be L/D)

% if just starting the recursion,
% initialize arrays to 0 
if (~exist('radii','var'))    
    radii = zeros(1,NT.nedge);        
end
if (isempty(radii))
    % initialize all arrays to 0
    radii = zeros(1,NT.nedge);        
end
if (~exist('stV','var'))    
    stV = zeros(1,NT.nedge);        
end
if (isempty(stV))
    % initialize all arrays to 0
    stV = zeros(1,NT.nedge);        
end



% set radius for current trunk
radii(trunkedge) = rtrunk;
ell0 = NT.edgelens(trunkedge);

%% junction node below this edge
junc = NT.edgenodes(trunkedge,2);

if (NT.degrees(junc)==1)
    % no daughter branches to deal with
    stV(trunkedge) = rtrunk^2*ell0;
    return

elseif (NT.degrees(junc)==2)
    error('not set up to deal with deg 2 nodes yet')
    
elseif (NT.degrees(junc)==3)
   
    edge1 = NT.nodeedges(junc,2);
    edge2 = NT.nodeedges(junc,3);
    
    % set radii for daughter branches    
    split21 = splitvals(edge2)/splitvals(edge1);
    r1 = ((rtrunk^a + rm^a)/(1+split21^(a/2)))^(1/a);
    r2 = r1*sqrt(split21);
    
    [radii,stV] = setRadiiWithRm(NT,edge1,a,rm,r1,splitvals,radii,stV);
    [radii,stV] = setRadiiWithRm(NT,edge2,a,rm,r2,splitvals,radii,stV);    
    
    % set volume of subtree with this trunk
    stV(trunkedge) = rtrunk^2*ell0 + stV(edge1)+stV(edge2);
else
    error('not set up to deal with node degrees > 3')
end



end