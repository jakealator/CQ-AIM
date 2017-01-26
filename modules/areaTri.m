function area = areaTri(x,y)

% From the 'shoelace' formula
area = 1/2*abs((x(1)-x(3))*(y(2)-y(1))-(x(1)-x(2))*(y(3)-y(1)));

end