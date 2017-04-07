
function absError = calculateAbsoluteError(M,dt,triAreas, qc, usSeries,us)

spaceErr = zeros(M,1);

for j=1:M
%     spaceErr(j) = norm(us(:,j)-usSeries(:,j),2)^2;
    spaceErr(j) = sum(triAreas.*abs(us(:,j)-usSeries(:,j)).^2);
end
absError = sqrt(dt*sum(spaceErr));

end
