% Linear interpolation function
function y = linearInterpolate(x0, y0, x1, y1, x)
    y = y0 + (y1 - y0) * ((x - x0) / (x1 - x0));
end
