%% sysprop
fid1 = fopen(strcat('inp_sysprop.dat'),'w');

% Write Mrb to file
line_str = num2str(size(Mrb,1)); % Number of lines to be used ==> determines number of DOF
fprintf(fid1,[line_str '\n']); 
for k1 = 1:size(Mrb,1)
    line_str = [num2str(Mrb(k1,1)) ' ' num2str(Mrb(k1,2))];
    fprintf(fid1,[line_str '\n']);
end

% Write Ainf to file
for k1 = 1:size(Add,1)
    line_str = [num2str(Add(k1,1)) ' ' num2str(Add(k1,2))];
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

%% simpar
fid2 = fopen(strcat('inp_simpar.dat'),'w');

% Write time parameters
line_str = num2str(ti);
fprintf(fid2,[line_str '\n']);
line_str = num2str(tf);
fprintf(fid2,[line_str '\n']);
line_str = num2str(dt);
fprintf(fid2,[line_str '\n']);
line_str = num2str(t_clutch);
fprintf(fid2,[line_str '\n']);

% Write initial states
% Initial position (eta0)
for k1 = 1:size(eta0,1)
    line_str = num2str(eta0(k1,1));
    fprintf(fid2,[line_str '\n']);
end

% Initial position (nu0)
for k1 = 1:size(nu0,1)
    line_str = num2str(nu0(k1,1));
    fprintf(fid2,[line_str '\n']);
end

% Write wind velocity
line_str = num2str(Uwind);
fprintf(fid2,[line_str '\n']);

% Write initial rotor speed
line_str = num2str(Omg_rt0);
fprintf(fid2,[line_str '\n']);

fclose(fid2);

%% rotorprop
fid3 = fopen(strcat('inp_rotorprop.dat'),'w');

% Write foils data to file
line_str = num2str(size(foilfile,1)); % Number of lines to be used ==> number of foils describing the blades
fprintf(fid3,[line_str '\n']); 
for k1 = 1:size(foilfile,1)
    line_str = char(foilfile(k1,1));
    fprintf(fid3,[line_str '\n']);
end

% Write number of blades, initial blade pitch angle, rotor radius, hub radius and hub height to file
line_str = num2str(Bl);
fprintf(fid3,[line_str '\n']); 
line_str = num2str(beta0);
fprintf(fid3,[line_str '\n']); 
line_str = num2str(Rtip);
fprintf(fid3,[line_str '\n']); 
line_str = num2str(Rhub);
fprintf(fid3,[line_str '\n']); 
line_str = num2str(Zhub);
fprintf(fid3,[line_str '\n']); 

% Write blade dimensions to file
line_str = num2str(size(Blade_dim,1)); % Number of lines to be used ==> number of elements describing the blades
fprintf(fid3,[line_str '\n']); 
for k1 = 1:size(Blade_dim,1)
    line_str = [num2str(Blade_dim(k1,1)) ' ' num2str(Blade_dim(k1,2)) ' '...
        num2str(Blade_dim(k1,3)) ' ' num2str(Blade_dim(k1,4)) ' ' num2str(Blade_dim(k1,5))];
    fprintf(fid3,[line_str '\n']);
end

fclose(fid3);

%% generator
fid4 = fopen(strcat('inp_generprop.dat'),'w');

% Write drivetrain properties to file
line_str = num2str(Igen);
fprintf(fid4,[line_str '\n']);
line_str = num2str(Ngr);
fprintf(fid4,[line_str '\n']);

% Write blade pitch controller properties to file
line_str = num2str(PC_KI);
fprintf(fid4,[line_str '\n']);
line_str = num2str(PC_KK);
fprintf(fid4,[line_str '\n']);
line_str = num2str(PC_KP);
fprintf(fid4,[line_str '\n']);
line_str = num2str(PC_MaxPit);
fprintf(fid4,[line_str '\n']);
line_str = num2str(PC_MaxRat);
fprintf(fid4,[line_str '\n']);
line_str = num2str(PC_MinPit);
fprintf(fid4,[line_str '\n']);
line_str = num2str(PC_RefSpd);
fprintf(fid4,[line_str '\n']);

% Write generator torque controller proerties to file
line_str = num2str(VS_CtInSp);
fprintf(fid4,[line_str '\n']);
line_str = num2str(VS_DT);
fprintf(fid4,[line_str '\n']);
line_str = num2str(VS_MaxRat);
fprintf(fid4,[line_str '\n']);
line_str = num2str(VS_MaxTq);
fprintf(fid4,[line_str '\n']);
line_str = num2str(VS_Rgn2K);
fprintf(fid4,[line_str '\n']);
line_str = num2str(VS_Rgn2Sp);
fprintf(fid4,[line_str '\n']);
line_str = num2str(VS_Rgn3MP);
fprintf(fid4,[line_str '\n']);
line_str = num2str(VS_RtGnSp);
fprintf(fid4,[line_str '\n']);
line_str = num2str(VS_RtPwr);
fprintf(fid4,[line_str '\n']);
line_str = num2str(VS_SlPc);
fprintf(fid4,[line_str '\n']);

fclose(fid4);

