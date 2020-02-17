function mask = makeMask(img)
% 10/09/2019
% Michael G. Moore, Michigan State University

% create a mask from a mean-image of fluoresence
% based on 10 percentile of intensity

img = (img - min(img(:)))/(max(img(:))-min(img(:)));
mask = img > .1;
mask = bwareaopen(mask,2000);
mask = imfill(mask,'holes');
mask = bwmorph(mask,'open');

end

