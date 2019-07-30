function strobeSeq = estimateStrobeSeq(imstack)

% strobeSeq = estimateStrobeSeq(imstack);
% estimates order of blue LED (1) and violet LED (2) by assuming that blue
% excitation images have higher RMS contrast. Has to be true for image
% stack overall

tic
fprintf('estimating LED strobe sequence...')
[nX,nY,nZ] = size(imstack);
imstack1   = imstack(ceil(nX*.2):floor(nX*.8),ceil(nY*.1):floor(nY*.8),1:2:nZ);
imstack2   = imstack(ceil(nX*.2):floor(nX*.8),ceil(nY*.1):floor(nY*.8),2:2:nZ);
seq1       = conditionDffMat(double(imstack1))';
seq2       = conditionDffMat(double(imstack2))';
seq1       = seq1./max(seq1(:));
seq2       = seq2./max(seq2(:));

m1         = repmat(mean(seq1),[size(seq1,1) 1]);
m2         = repmat(mean(seq2),[size(seq2,1) 1]);
rms1       = mean(sqrt(sum((seq1-m1).^2)/size(seq1,1)));
rms2       = mean(sqrt(sum((seq2-m2).^2)/size(seq2,1)));

if rms1 > rms2; strobeSeq = [1 2]; else strobeSeq = [2 1]; end
fprintf(' done after %1.1f min\n',toc/60)