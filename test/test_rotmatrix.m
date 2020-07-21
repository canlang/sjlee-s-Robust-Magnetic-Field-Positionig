R1 = arrayfun(@(x)(rotMat(:,:,1)*[cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]),...
        ps.mag_heading,'UniformOutput',false); 

R2 = arrayfun(@(x) euler2rotMat(euler(1,1),euler(1,2),x),...
    ps.mag_heading,'UniformOutput',false);
    
    
R3 = arrayfun(@(x) euler2rotMat(euler(1,1),euler(1,2),x),...
        -euler(1,3)-ps.mag_heading,'UniformOutput',false);
rotatedMag3 = cell2mat(cellfun(@(x)(x*tM(i,:)')',R,'UniformOutput',false));

R4 = arrayfun(@(x) euler2rotMat(euler(1,1),euler(1,2),x),...
        euler(1,3)+ps.mag_heading,'UniformOutput',false);
rotatedMag4 = cell2mat(cellfun(@(x)(x.'*tM(i,:)')',R,'UniformOutput',false));
    
