function plotTimeSolution(meshStructExt,xExt,uExt,meshStruct,x,u,M,movieFlag) 

if movieFlag
    figure
    for j=1:M
        pltsln(meshStructExt,xExt,uExt(:,j))
        view(2)
        colormap(jet)
        axis([-1 1 -1 1 -1.1 1.8 -1.1 1.8])
        axis equal
        hold on
        pltsln(meshStruct,x,u(:,j))
        view(2)
        colormap(jet)
        axis equal
        title(j)
        pause(.1)
        hold off
    end
else
    figure
    j=1;
    pltsln(meshStructExt,xExt,uExt(:,j))
    view(2)
    colormap(jet)
    axis([-1 1 -1 1 -1.1 1.8 -1.1 1.8])
    axis equal
    hold on
    pltsln(meshStruct,x,u(:,j))
    view(2)
    colormap(jet)
    axis equal
    title('Time t=0')
    
    figure
    j=18;
    pltsln(meshStructExt,xExt,uExt(:,j))
    view(2)
    colormap(jet)
    axis([-1 1 -1 1 -1.1 1.8 -1.1 1.8])
    axis equal
    hold on
    pltsln(meshStruct,x,u(:,j))
    view(2)
    colormap(jet)
    axis equal
    title('Time t=1.15')
    
    figure
    j=36;
    pltsln(meshStructExt,xExt,uExt(:,j))
    view(2)
    colormap(jet)
    axis([-1 1 -1 1 -1.1 1.8 -1.1 1.8])
    axis equal
    hold on
    pltsln(meshStruct,x,u(:,j))
    view(2)
    colormap(jet)
    axis equal
    title('Time t=2.3')
    
    
    figure
    j=54;
    pltsln(meshStructExt,xExt,uExt(:,j))
    view(2)
    colormap(jet)
    axis([-1 1 -1 1 -1.1 1.8 -1.1 1.8])
    axis equal
    hold on
    pltsln(meshStruct,x,u(:,j))
    view(2)
    colormap(jet)
    axis equal
    title('Time t=3.5')
    
    figure
    j=70;
    pltsln(meshStructExt,xExt,uExt(:,j))
    view(2)
    colormap(jet)
    axis([-1 1 -1 1 -1.1 1.8 -1.1 1.8])
    axis equal
    hold on
    pltsln(meshStruct,x,u(:,j))
    view(2)
    colormap(jet)
    axis equal
    title('Time t=4.5')
    
end

end
