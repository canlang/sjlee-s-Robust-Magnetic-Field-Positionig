R1 = arrayfun(@(x)(rotMat(:,:,1)*[cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]),...
        ps.mag_heading,'UniformOutput',false); 

R2 = arrayfun(@(x) euler2rotMat(euler(1,1),euler(1,2),x),...
    ps.mag_heading,'UniformOutput',false);

R3 = arrayfun(@(x) [cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]*rotMat(:,:,i)...
        ,-ps.mag_heading,'UniformOutput',false); 
    
    
rotatedMag1 = cell2mat(cellfun(@(x)(x*tM(1,:)')',R,'UniformOutput',false));


rotatedMag2 = cell2mat(cellfun(@(x)(x.'*tM(1,:)')',R,'UniformOutput',false));    

