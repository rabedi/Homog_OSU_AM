%fid = fopen('grain1.inp', 'r');
fid = fopen('grain1_complete.inp', 'r');
g1 = Grain;
g1 = g1.read(fid, 1);
save('g1Output.mat', 'g1');
a = [g1.grainxm, g1.grainxM, g1.grainym, g1.grainyM];
fclose(fid);
%listDomainEdges_grain1 = g1.getDomainEdges();

fido = fopen('grain1_out.inp', 'w');
g1.write(fido);
fclose(fido);

g2 = Grain;
fidi = fopen('grain1_out.inp', 'r');
g2 = g2.read(fidi, 1);
figure(1);
%g2.plotMesh();
g2.fillMesh(); %clr_tri, clr_quad)
xlim([260,272]);
ylim([5,10]);
axis('square');


% grain is inside x [2, 272]  times y [1, 161.47]
%xmin = 0, xmax = 10, ymin = 0, ymax = 10;
xmin = 260; xmax = 280; ymin = 6; ymax = 8;
renumbering = 1;
[g1Intersect, hasIntersection] = g1.Extract(xmin, xmax, ymin, ymax, renumbering);
fido = fopen('grain1Intersect_out.inp', 'w');
g1Intersect.write(fido);
fclose(fido);
figure(2);
g1Intersect.plotMesh();
listDomainEdges_grain1_intersected = g1Intersect.getDomainEdges();
a = 12;
