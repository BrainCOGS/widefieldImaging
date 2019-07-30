function comVec = COM(mat,taxis)

% comMat = COM(mat,taxis)
% calculates center of mass for vectors going along rows of mat

if taxis(1) == 0
    taxis = taxis + mode(diff(taxis));
end

% add offset as negative numbers mess up COM
mat    = mat+repmat(abs(min(mat,[],2)),[1 size(mat,2)]);    

comVec = ((taxis*mat')./sum(mat'))';
