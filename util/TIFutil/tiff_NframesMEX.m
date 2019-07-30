function nZ = tiff_NframesMEX(filename)

temp = ecs.imfinfox(filename);
nZ   = temp.fileFrames;