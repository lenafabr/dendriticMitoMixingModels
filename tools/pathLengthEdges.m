function lend = pathLengthEdges(NST)
% compute the path length to the distal tip of each branch in a tree


lend=zeros(NST.nedge,1);
for ec=1:NST.nedge
    lend(ec)=NST.edgelens(ec,1);
    n1=NST.edgenodes(ec,1);
    degp=NST.degrees(n1,1);
    while(degp ~= 1)
        ecp=NST.nodeedges(n1,1);
  	  lend(ec)=lend(ec) +NST.edgelens(ecp) ;
      n1=NST.edgenodes(ecp,1);
      degp=NST.degrees(n1,1);
    end
end

end