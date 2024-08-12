fid = fopen('grain1 (1).inp', 'r');
g1 = Grain;
g1 = g1.read(fid, 1)
a = [g1.grainxm, g1.grainxM, g1.grainym, g1.grainyM]
fclose(fid);

fido = fopen('grain_this_the_one(9).inp', 'w');
g1.write(fido);
fclose(fido);
g1.plotMesh


% grain is inside x [2, 272]  times y [1, 161.47]
xmin = 0; xmax = 10; ymin = 0; ymax = 10;
[g1Intersect, hasIntersection] = g1.Extract(xmin, xmax, ymin, ymax);
fido = fopen('grain1Intersect_out(test5).inp', 'w');
g1Intersect.write(fido);
fclose(fido);
g1Intersect.plotMeshRestricted


