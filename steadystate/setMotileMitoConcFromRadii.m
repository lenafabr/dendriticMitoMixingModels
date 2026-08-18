function [rhoWvals] = setMotileMitoConcFromRadii(NT,trunkedge,rhoWtrunk,radii,rhoWvals)
% walk down a tree, starting with the trunkedge
% using precalculated radii 
% set the linear concentrations of motile mitochondria (rhoWvals)
% ** assuming walking mitos split in proportion to r^2 **
% NT = network object
% trunkedge = index for the trunk edge
% rhoWtrunk = motile density in the trunk
% radii = radius for each edge

% if just starting the recursion,
% initialize arrays to 0 
if (~exist('rhoWvals','var'))    
    rhoWvals = zeros(1,NT.nedge);        
end
if (isempty(rhoWvals))
    % initialize all arrays to 0
    rhoWvals = zeros(1,NT.nedge);        
end


% set rho for current trunk
rhoWvals(trunkedge) = rhoWtrunk;

%% junction node below this edge
junc = NT.edgenodes(trunkedge,2);

if (NT.degrees(junc)==1)
    % no daughter branches to deal with
    % done with the recurrsion
    return

elseif (NT.degrees(junc)==2)
    % same linear density just below this node
    edge1 = NT.nodeedges(junc,2);
    rhoWvals = setMotileMitoConcFromRadii(NT,edge1,rhoWvals(trunkedge),radii,rhoWvals);    
elseif (NT.degrees(junc)==3)
   
    edge1 = NT.nodeedges(junc,2);
    edge2 = NT.nodeedges(junc,3);
    
    % get linear densities for daughter branches
    r1 = radii(edge1); r2 = radii(edge2);
    rhoW1 = r1^2/(r1^2+r2^2)*rhoWvals(trunkedge);
    rhoW2 = r2^2/(r1^2+r2^2)*rhoWvals(trunkedge);
    
    [rhoWvals] = setMotileMitoConcFromRadii(NT,edge1,rhoW1,radii,rhoWvals);
    [rhoWvals] = setMotileMitoConcFromRadii(NT,edge2,rhoW2,radii,rhoWvals);        
    
else
    error('not set up to deal with node degrees > 3')
end



end