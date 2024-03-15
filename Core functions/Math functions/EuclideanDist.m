function [dist] = EuclideanDist(V1,V2)
%EUCLEANDIST Summary of this function goes here
%computes Euclidean distance between two vectors
%   Detailed explanation goes here
dist=sqrt(sum((V1-V2).^2));
end

