
function relativeError = calculateRelativeError(M,usSeries,us)

spaceNormExact = zeros(M,1);
spaceErr = zeros(M,1);

for j=1:M
    spaceNormExact(j) = norm(abs(usSeries(:,j)),2)^2;
    spaceErr(j) = norm(abs(us(:,j)-usSeries(:,j)),2)^2;
end
spaceTimeErr = sum(spaceErr);
spaceTimeNorm = sum(spaceNormExact);

relativeError = sqrt(spaceTimeErr/spaceTimeNorm);

end