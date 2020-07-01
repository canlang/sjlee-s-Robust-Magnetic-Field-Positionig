function [distance, index] = findNearestLocation(locs_1, locs_2)
[distance,index] = pdist2(locs_1,locs_2,'euclidean','Smallest',1);