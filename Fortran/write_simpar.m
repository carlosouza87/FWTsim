fid2 = fopen(strcat('inp_simpar.txt'),'w');

% Write time parameters
line_str = num2str(ti);
fprintf(fid1,[line_str '\n']);
line_str = num2str(tf);
fprintf(fid1,[line_str '\n']);
line_str = num2str(dt);
fprintf(fid1,[line_str '\n']);

% Write initial states
line_str = num2str(size(eta0,1)); % Number of lines to be used ==> determines number of DOF
fprintf(fid1,[line_str '\n']); 

% Initial position (eta0)
for k1 = 1:size(eta0,1)
    line_str = num2str(eta0(k1,1));
    fprintf(fid1,[line_str '\n']);
end

% Initial position (nu0)
for k1 = 1:size(nu0,1)
    line_str = num2str(nu0(k1,1));
    fprintf(fid1,[line_str '\n']);
end

fclose(fid2);
