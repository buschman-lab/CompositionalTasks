function [ClusterLabels,Linkage,squarFormY] = HierarchicalClustring(input,MaxClust)
%HIERARCHICALCLUSTRING Summary of this function goes here
%   Detailed explanation goes here

Y=pdist(input);
squarFormY=squareform(Y);
Linkage=linkage(Y,'ward');
ClusterLabels = cluster(Linkage,'maxclust',MaxClust);

end

