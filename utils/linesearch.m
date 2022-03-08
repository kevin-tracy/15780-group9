
function alpha = linesearch(x,dx)
    alpha = 1.0;
    for i = 1:length(x)
        if dx(i)<0
            alpha = min(alpha,-x(i)/dx(i));
        end
    end
end