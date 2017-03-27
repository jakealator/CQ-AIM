function y=FX(w,a,b,c)
    y=sqrt(pi/b)*sin(a*((1i*w)/(2*b)+c))*...
        exp(-1/(4*b)*(a^2-4*1i*b*c*w+w.^2));
end