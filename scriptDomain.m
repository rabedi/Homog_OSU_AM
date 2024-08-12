dom = Domain;
fid = fopen('domain.mat', 'r');
if (fid > 0)
    fclose(fid);
    dom = load('domain.mat');
    dom = dom.dom;
else
    fid = fopen('Final_mesh.inp', 'r');
    tic
    dom = dom.readBeforeMaterialAssembly(fid);
    toc
    save('domain.mat', 'dom');
    fclose(fid);
end

[dom.domainxm, dom.domainxM, dom.domainym, dom.domainyM] = ...
dom.ComputeDomainBounds();
dom.grains{1}.fillMesh();
save('grain1Alone.fig');
print('-dpng', 'grain1Alone.png');

figure(2);
dom.grains{223}.fillMesh();
save('grain223Alone.fig');
print('-dpng', 'grain223Alone.png');
xlim([3270, 3310]);
ylim([2960,3007]);
save('grain223AloneZoom.fig');
print('-dpng', 'grain223AloneZoom.png');


domIntersect = Domain;
%xmin = 260; xmax = 280; ymin = 6; ymax = 8;
xmin = 0; xmax = 500; ymin = 0; ymax = 200;
renumbering = 1;
tol = 1e-5;
domIntersect = dom.Extract(xmin, xmax, ymin, ymax, renumbering, tol);
[vecSegmentXss, vecSegmentYss, vecSegmentNodeIDss, grainIDs] = domIntersect.PlotGetDomainEdges();
fileNameBaseWithSVENo = 'intersectedSVE';
configuration = OSU_SpecificOutputConfiguration;
center_4DispCal = [];
domIntersect.writeDifferentBC_SVEs(fileNameBaseWithSVENo, configuration, center_4DispCal);

a = 12;

% %dom.fill_visualize(1);
% %save('grain1p.fig');
% %print('-dpng', 'grain1p.png');
% tic
% tol = 1e-6;
% dom = dom.FormAllGrainEdgeNeighborhoodInfo(tol);
% %vecGrainBoundaries = dom.getVecGrainBoundaries();
% toc
% a = 12;
