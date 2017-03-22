function uiHat = generateUIHat(uiHatFun, femStruct)

    % Define incident field on centroids
     uiHat = uiHatFun(femStruct.centroids);

end