function writeTiffStack(imstack,fn)

% write tiff stack from 3d matlab matrix

tic
fprintf('writing tiff stack...')
if ~isinteger(imstack)
    imstack = uint16(imstack);
end

imwrite(imstack(:,:,1),fn)
for k = 2:size(imstack,3)
    imwrite(imstack(:,:,k),fn, 'writemode', 'append');
    if rem(k,500) == 0; fprintf('\n   %d/%d frames written',k,size(imstack,3)); end
end
t = toc;
fprintf('\n   done after %1.2f minutes\n',t/60)