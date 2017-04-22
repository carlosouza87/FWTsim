fid1 = fopen(strcat('inp_rotor.dat'),'w');

% Write foils data to file
line_str = num2str(size(foilfile,1)); % Number of lines to be used ==> number of foils describing the blades
fprintf(fid1,[line_str '\n']); 
for k1 = 1:size(foilfile,1)
    line_str = char(foilfile(k1,1));
    fprintf(fid1,[line_str '\n']);
end

% Write number of blades, rotor radius, hub radius and hub height to file
line_str = num2str(Br);
fprintf(fid1,[line_str '\n']); 
line_str = num2str(Rtip);
fprintf(fid1,[line_str '\n']); 
line_str = num2str(Rhub);
fprintf(fid1,[line_str '\n']); 
line_str = num2str(Zhub);
fprintf(fid1,[line_str '\n']); 

% Write blade dimensions to file
line_str = num2str(size(Blade_dim,1)); % Number of lines to be used ==> number of elements describing the blades
fprintf(fid1,[line_str '\n']); 
for k1 = 1:size(Blade_dim,1)
    line_str = [num2str(Blade_dim(k1,1)) ' ' num2str(Blade_dim(k1,2)) ' '...
        num2str(Blade_dim(k1,3)) ' ' num2str(Blade_dim(k1,4)) ' ' num2str(Blade_dim(k1,5))];
    fprintf(fid1,[line_str '\n']);
end

fclose(fid1);

