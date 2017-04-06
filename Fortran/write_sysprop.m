fid1 = fopen(strcat('inp_sysprop.txt'),'w');

% Write Mrb to file
for k1 = 1:size(Mrb,1)
    line_str = [num2str(Mrb(k1,1)) ' ' num2str(Mrb(k1,2))];
    fprintf(fid1,[line_str '\n']);
end

% Write Ainf to file
for k1 = 1:size(Ainf,1)
    line_str = [num2str(Ainf(k1,1)) ' ' num2str(Ainf(k1,2))];
    fprintf(fid1,[line_str '\n']);
end

% Write Dl to file
for k1 = 1:size(Dl,1)
    line_str = [num2str(Dl(k1,1)) ' ' num2str(Dl(k1,2))];
    fprintf(fid1,[line_str '\n']);
end

% Write Dq to file
for k1 = 1:size(Dq,1)
    line_str = [num2str(Dq(k1,1)) ' ' num2str(Dq(k1,2))];
    fprintf(fid1,[line_str '\n']);
end

% Write Chs to file
for k1 = 1:size(Chs,1)
    line_str = [num2str(Chs(k1,1)) ' ' num2str(Chs(k1,2))];
    fprintf(fid1,[line_str '\n']);
end

% Write Cmr to file
for k1 = 1:size(Cmr,1)
    line_str = [num2str(Cmr(k1,1)) ' ' num2str(Cmr(k1,2))];
    fprintf(fid1,[line_str '\n']);
end

fclose(fid1);