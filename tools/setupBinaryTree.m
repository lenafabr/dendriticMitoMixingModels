function NT=setupAxonTree(nLevels,l,theta_list,scale)
% Build a binary tree
% nLevels = number of branching levels
% l = branch length
% theta_list = angles of branches relative to parent (suggest [pi/4, pi/6])
% scale = scaling factor for branch length going down the tree?

%% Initialize root node

nodepos = [-l,0;0 0];  % root horizontal line
edgenodes = [1 2];     % zeroth edge
node_counter = 2;

% Each branch: [start_x, start_y, direction_angle, length, parent_index]
branches = [0, 0, 0, l, 2];  % start from root edge (angle=0 along x-axis)

%% Build tree
for level = 1:nLevels
    new_branches = [];
    if level==1
            theta=theta_list(1);
    else
            theta=theta_list(2);
    end
    for i = 1:size(branches,1)
        x = branches(i,1);
        y = branches(i,2);
        ang = branches(i,3);
        len = branches(i,4);
        parent = branches(i,5);
        

        % Left child (relative angle +theta)
        angL = ang + theta;
        xL = x + len*cos(angL);
        yL = y + len*sin(angL);
        nodepos(end+1,:) = [xL, yL];
        idxL = size(nodepos,1);
        edgenodes(end+1,:) = [parent, idxL];
        new_branches = [new_branches; xL, yL, angL, len*scale, idxL];
        
        % Right child (relative angle -theta)
        angR = ang - theta;
        xR = x + len*cos(angR);
        yR = y + len*sin(angR);
        nodepos(end+1,:) = [xR, yR];
        idxR = size(nodepos,1);
        edgenodes(end+1,:) = [parent, idxR];
        new_branches = [new_branches; xR, yR, angR, len*scale, idxR];
    end
    
    branches = new_branches;  % next level    
end

% convert to network object
NT = NetworkObj();
NT.nodepos = nodepos;
NT.edgenodes = edgenodes;
NT.setupNetwork();

% make network object a directed tree
NT.rootnode = 1;
directedTreeEdges(NT,NT.rootnode,false(1,NT.nedge),false(1,NT.nedge));
