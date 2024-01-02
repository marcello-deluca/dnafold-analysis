function [points] = MakeBeautifiedFrame(FRAME, RADIUS, COLOR)
    spline = cscvn(transpose(FRAME));
    [points] = fnplt(spline,'r',2);
    spline2 = cscvn(points);
    [points] = fnplt(spline2, 'r',2);
    spline3 = cscvn(points);
    [points] = fnplt(spline3, 'r',2);
    tubeplot(points,RADIUS,128,[],COLOR);
    %camlight;
    view(-62,5);
   % camlight;
    view(-45,5);
    %camlight;
    view(-62,5);
    %light
    set(gcf, 'color', 'white');
end