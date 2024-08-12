function CopyFileContent(fidwrite, fileNameRead)
fidread = fopen(fileNameRead, 'r');
if (fidread < 0)
    fprintf(1, 'Cannot open %s\n', fileNameRead);
    pause;
end
tline = fgetl(fidread);
while ~feof(fidread)
    fprintf(fidwrite, '%s\n', tline);
    tline = fgetl(fidread);
end
if (length(tline) > 0)
    fprintf(fidwrite, '%s\n', tline);
end    
fclose(fidread);