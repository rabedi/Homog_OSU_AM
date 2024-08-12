classdef OSU_SpecificOutputConfiguration
    properties
%        E = 1000;
        nu = 0.3;
%        density = 1e-5;
        isElastic = 1;
%        epsMagnitude = 1e-4;

        E = 175;
        density = 7000 * 1e-12;
        epsMagnitude = 0.02;
    end
    methods
        function print_segment(obj, fidwrite, fileNameReadWOExt)
            fileNameRead = [fileNameReadWOExt, '.txt']; 
            CopyFileContent(fidwrite, fileNameRead);
        end
        function print_grain_material_properties(obj, grainID, fidwrite)
            fprintf(fidwrite, '*MATERIAL, NAME=MATERIAL_%d\n', grainID);
            fprintf(fidwrite, '*DENSITY\n');
            fprintf(fidwrite, '\t%d\n', obj.density);
            fprintf(fidwrite, '**\n');
            fprintf(fidwrite, '*ELASTIC, TYPE=ISO\n');
            fprintf(fidwrite, '\t%d,\t%d\n', obj.E, obj.nu);
            fprintf(fidwrite, '**\n');
        end
        
        function [eps_xx_eps_yy_eps_xyS, addedNames] = getStrains(obj)
                eps_xx_eps_yy_eps_xyS = obj.epsMagnitude * [1, 0, 0; 0, 1, 0; 0, 0, 0.5];
                addedNames = {'_1', '_2', '_3'};
        end
    end
end
