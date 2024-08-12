%scale in x, y
scales = [0.2862    0.2581];
% elastic modulus
E = 175; % horizontal is 172 GPa, vertical 60 GPa, we start with this ...
rho = 7000 * 1e-12; % factor of 1e-12 was determined to map to target rho
strengthScale = 1076e-3; % vertical yield strength in MPa -> GPa needs factor 1e-3
strainScale = strengthScale / E;
% to Ryan: If the loading is small (not causing softening and unloading) change 3 ->
% larger values
maxAppliedStrain = 3.0 * strainScale; % we apply strength about 3x that what is needed to get to yield strength
%0.02 is used for strain magnitude;

%%%%% cohesive model parameters
% there are about 20 printed layer -> spacing about cohesive surfaces is
% spacing = height / 20 = 851/20 =~ 42

%% initial stiffness of cohesive model K, spacing -> efective compliance = 1/K.spacing
% we want the effective compliance ~ 160 in vertical direction
% .. Compliance_vertical_effective = 1/E + 1/K spacing ->
% 1/160 = 1/175 + 1/k .spacing -> K . spacing = 1867, spacing = 42 -> K =
% 45

% to Ryan: change K here
coh_K = 32;
% to Ryan: sigma_max for cohesive model
coh_strength = 1080e-3; % 1.08

%computed from this coh_energy = coh_strength * base, assuming right
%segment (unloading) is say r2l times left segment, energy = coh_strengh /
%K / 2 * (1 + r2l)
r2l = 2;
coh_energy = coh_strength / coh_K / 2.0 * (1 + r2l)
% use
% to Ryan: If you want, you can put your value of coh_energy here
coh_energy = 0.046;


%%%%%%%%%%%%%%%%%%%%%%%%



% the domain is [2, 3303] x [1, 3000], grid points are locations 10i, 10j
% (initial grid at spacing 10) - number below specifies how larger SVEs are
% relative to 10
% to Ryan: you can increase the SVE size from 10 to something larger say 100
% but processing of SVEs take longer ...
SVE_size_multipleOf10 = 100; % 100 see if it works
% if 1: initial band size 8 in x, 9 in y is included in forming SVEs. 
% if 0: domain starts from [10, 10]
includeSmallerSVEsLeftBottom = 1; 

renumbering = 1;
tol = 1e-4;
plotSVE_Bndry = 1;


configuration = OSU_SpecificOutputConfiguration;
configuration.epsMagnitude = 1e-4;

center_4DispCal = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
runSign = ['OSU_sm_', num2str(includeSmallerSVEsLeftBottom), '_size_', num2str(SVE_size_multipleOf10)];
fidSum = fopen([runSign, '_summary.txt'], 'w');

SVE_size = 10 * SVE_size_multipleOf10;
if (includeSmallerSVEsLeftBottom)
    SVE_xlimts = 0:SVE_size:3303;
    SVE_ylimts = 0:SVE_size:3300;
else
    SVE_xlimts = 10:SVE_size:3303;
    SVE_ylimts = 10:SVE_size:3300;
%    SVE_xlimts = SVE_size:SVE_size:3303;
%    SVE_ylimts = SVE_size:SVE_size:3300;
end
num_SVEs_x = length(SVE_xlimts) - 1;
num_SVEs_y = length(SVE_ylimts) - 1;
num_SVEs = num_SVEs_x * num_SVEs_y; 
fprintf(fidSum, 'num_SVEs_x\t%d\tnum_SVEs_y\t%d\tnum_SVEs\t%d\n', num_SVEs_x, num_SVEs_y, num_SVEs);
xmins = zeros(num_SVEs_x, num_SVEs_y);
xmaxs = zeros(num_SVEs_x, num_SVEs_y);
ymins = zeros(num_SVEs_x, num_SVEs_y);
ymaxs = zeros(num_SVEs_x, num_SVEs_y);
sveNums = zeros(num_SVEs, 2);
sveNames = cell(num_SVEs, 1);
cntr = 0;
for i = 1:num_SVEs_x
    xm = SVE_xlimts(i);
    xM = SVE_xlimts(i + 1);
    for j = 1:num_SVEs_y
        ym = SVE_ylimts(j);
        yM = SVE_ylimts(j + 1);
        xmins(i, j) = xm;
        xmaxs(i, j) = xM;
        ymins(i, j) = ym;
        ymaxs(i, j) = yM;
        cntr = cntr + 1;
        sveNums(cntr, 1) = i;
        sveNums(cntr, 2) = j;
        name = ['SVE_', num2str(cntr)];
        SVEName = [runSign, '_', name];
        sveNames{cntr} = SVEName;
        fprintf(fidSum, 'cntr\t%d\ti\t%d\tj\t%d\txm\t%f\txM\t%f\tym\t%f\tyM\t%f\tSVEName\t%s\n', cntr, i,j, xm, xM, ym, yM, SVEName);
    end
end
fclose(fidSum);


%%%%%%%%% reading the domain or recreating it, as needed
dom = Domain;
fid = fopen('domain.mat', 'r');
if (fid > 0)
    fclose(fid);
    dom = load('domain.mat');
    dom = dom.dom;
else
    % 
    fid = fopen('Final_mesh.inp', 'r');
    tic
    dom = dom.readBeforeMaterialAssembly(fid);
    toc
    save('domain.mat', 'dom');
    fclose(fid);
end

%%%%% print scaled domain
if 0
    % first write scaled dom, ...
    dom_scaled = dom;
    fid = fopen('Final_mesh_scaled.inp', 'w');
    configuration = OSU_SpecificOutputConfiguration;
    configuration.print_segment(fid, '__part_0_header');
    
    dom_scaled.domainxm = dom_scaled.domainxm * scales(1);
    dom_scaled.domainxM = dom_scaled.domainxM * scales(1);
    dom_scaled.domainym = dom_scaled.domainym * scales(2);
    dom_scaled.domainyM = dom_scaled.domainyM * scales(2);
    for grindi = 1:dom_scaled.num_grains
        num_node = dom_scaled.grains{grindi}.nodeNum;
        for j = 1:num_node
            tmp = dom_scaled.grains{grindi}.nodeCrds(j, 1);
            dom_scaled.grains{grindi}.nodeCrds(j, 1) = tmp * scales(1);
            tmp = dom_scaled.grains{grindi}.nodeCrds(j, 2);
            dom_scaled.grains{grindi}.nodeCrds(j, 2) = tmp * scales(2);
        end
        fprintf(fid, '*Part, name=GRAIN-%d\n', dom_scaled.grains{grindi}.id);
        dom_scaled.grains{grindi}.write(fid);
    end
    fclose(fid);
end

cntr = 0;

fidEmpyu = fopen([runSign, '_empty.txt'], 'w');

for i = 1:num_SVEs_x
    xm = SVE_xlimts(i);
    xM = SVE_xlimts(i + 1);
    for j = 1:num_SVEs_y
        cntr = cntr + 1;
        SVEName = sveNames{cntr};
        fprintf(1, '%d->SVE(%d,%d)\n', cntr, i, j); 

        domIntersect = Domain;
        xm = xmins(i, j);
        xM = xmaxs(i, j);
        ym = ymins(i, j);
        yM = ymaxs(i, j);
        domIntersect = dom.Extract(xm, xM, ym, yM, renumbering, tol, scales);
        if (domIntersect.num_grains == 0)
            fprintf(fidEmpyu, '%d\t%d\t%d\n', cntr, i, j);
            continue;
        end
        if (plotSVE_Bndry)
            h = figure(1);
            SVEName_plot = [SVEName, '_boundary'];
            [vecSegmentXss, vecSegmentYss, vecSegmentNodeIDss, grainIDs] = domIntersect.PlotGetDomainEdges();
            print('-dpng', [SVEName_plot, '.png']);
            %save([SVEName_plot, '.fig']);
            close(h);
        end
        domIntersect.writeDifferentBC_SVEs(SVEName, configuration, center_4DispCal);
        clear domIntersect;
    end
end
fclose(fidEmpyu);

% [dom.domainxm, dom.domainxM, dom.domainym, dom.domainyM] = ...
% dom.ComputeDomainBounds();
% dom.grains{1}.fillMesh();
% save('grain1Alone.fig');
% print('-dpng', 'grain1Alone.png');
% figure(2);
% dom.grains{223}.fillMesh();
% save('grain223Alone.fig');
% print('-dpng', 'grain223Alone.png');
% xlim([3270, 3310]);
% ylim([2960,3007]);
% save('grain223AloneZoom.fig');
% print('-dpng', 'grain223AloneZoom.png');
