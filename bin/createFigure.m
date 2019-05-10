function f = createFigure(left,bottom,width,height)
% create a figure with handle 'f' and specify size and position of figure
% units are given in fractions of the screen dimensions

sz = get(groot,'ScreenSize');


f = figure;
set(f,'OuterPosition',[left*sz(3) bottom*sz(4) width*sz(3) height*sz(4)])

end