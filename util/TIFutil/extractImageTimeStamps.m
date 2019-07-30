function ts = extractImageTimeStamps(tiffFile)

% reads tiff stack headers and returns time stamps in seconds
fprintf('extracting time stamps from tiff...')
if iscell(tiffFile)
    ts = [];
    for ii = 1:length(tiffFile)
        info = parseTIFFHeaderLP(tiffFile{ii});
        close
        tst  = zeros(length(info),1);    
        for jj = 1:length(tst)
            thist          = info(jj).Time_From_Start;
            [~,~,~,H,MN,S] = datevec(thist);
            tst(jj) = S + 60*MN + 3600*H;
        end
        ts = [ts; tst];
    end
else
    info = parseTIFFHeaderLP(tiffFile);
    close
    ts   = zeros(length(info),1);
    
    for ii = 1:length(ts)
        thist          = info(ii).Time_From_Start;
        [~,~,~,H,MN,S] = datevec(thist);
        ts(ii) = S + 60*MN + 3600*H;
    end
end
fprintf('\n')