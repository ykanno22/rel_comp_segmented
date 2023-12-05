clear;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some parameters
% --->
% param.kappa = 10^(-5);
param.kappa = 2.0 * 10^(0);
param.big_M = 10^(3);
%
% Flag.exe = 1; %%%% read sol-file and create figure
% Flag.exe = 2; %%%% create lp-file
Flag.exe = 3; %%%% solve problem
% Flag.miqcqstrat = 1;  %%% QP-base
Flag.miqcqstrat = 2;  %%% LP-base
% Flag.display = 0;
Flag.display = 1;

load('tri-modulus_1d_data_set.mat');
num.data = length(list_of_noisy_x);
num.dim  = 2;
num.b_point = 5;
nr = num.data;
nb = num.b_point;  %%% #break points
%
data_eps = list_of_noisy_x';
data_sig = list_of_noisy_f';
clear list_of_noisy_x list_of_noisy_f
%
filename.std = 'miqp_segmented_reg';
% <---
% Some parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discrete optimization
% --->
fprintf('\n  ##### MIQP >>>>> \n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename.lp = strcat(filename.std,'.lp');
filename.cplex_log = strcat(filename.std,'_cplex.log');
filename.cplex_sol = strcat(filename.std,'_cplex.sol');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Define variable types -->
%%%% variables : [v]      \in (nb*nr)                              %%%%
%%%% variables : [a]      \in (nb)                                 %%%%
%%%% variables : [b]      \in (nb)                                 %%%%
%%%% variables : [c]      \in (nb)                                 %%%%
%%%% variables : [t]      \in (nb*nr)                              %%%%
char.v = repmat('C',nb*nr,1);
char.a = repmat('C',nb,1);
char.b = repmat('C',nb,1);
char.c = repmat('C',nb,1);
char.t = repmat('B',nb*nr,1);
char.c = repmat('C',nb,1);
char.t = repmat('B',nb*nr,1);

vartype =...
    [char.v; char.a; char.b; char.c; char.t];
clear char

num.all_var = length(vartype);
%%%% <-- Define variable types
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective function
% --->
obj.v = repmat(1,nb*nr,1);
obj.a = repmat(0,nb,1);
obj.b = repmat(0,nb,1);
obj.c = repmat(0,nb,1);
obj.t = repmat(0,nb*nr,1);
for i=1:nb
    obj.t(i*nr) = param.kappa;
end

vec_obj =...
    sparse([obj.v; obj.a; obj.b; obj.c; obj.t]);
clear obj
% <---
% Objective function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define lower/upper bounds
% --->
%%%% lower bounds
lb.v = repmat(0,nb*nr,1);
lb.a = repmat(-Inf,nb,1);
lb.b = repmat(-Inf,nb,1);
lb.c = repmat(-Inf,nb,1);
lb.t = repmat(0,nb*nr,1);
%
vec_lb = [lb.v; lb.a; lb.b; lb.c; lb.t];
clear lb

%%%% upper bounds
ub.v = repmat(Inf,nb*nr,1);
ub.a = repmat(Inf,nb,1);
ub.b = repmat(Inf,nb,1);
ub.c = repmat(Inf,nb,1);
ub.t = repmat(1,nb*nr,1);
%
vec_ub = [ub.v; ub.a; ub.b; ub.c; ub.t];
clear ub
% <---
% Define lower/upper bounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear inequalities
% --->
ctype = [];
con_ax = [];
con_bx = [];

%%% inequalities for (t) : monotonicity (1)
pp = mk_pointer_1(num);
%
A_t = zeros(nb*(nr-1),num.all_var);
%
A_t_unit = eye(nr-1);
for l=1:(nr-1)
    A_t_unit(l,l+1) = -1;
end
qq = 0;
for i=1:nb
    A_t((qq+1):(qq+nr-1),(pp.t+1):(pp.t+nr)) = A_t_unit;
    qq = qq + nr - 1;
    pp.t = pp.t + nr;
end
A_t = sparse(A_t);
b_t = sparse(nb*(nr-1),1);
%
con_ax = [con_ax; A_t];
con_bx = [con_bx; b_t];
ctype  = [ctype, repmat('U',1,nb*(nr-1))];
clear pp qq A_t b_t

%%% inequalities for (t) : monotonicity (2)
pp = mk_pointer_1(num);
%
A_t = zeros((nb-1)*nr,num.all_var);
%
A_t_unit = [-eye(nr), eye(nr)];
qq = 0;
for i=1:(nb-1)
    A_t((qq+1):(qq+nr),(pp.t+1):(pp.t+(2*nr))) = A_t_unit;
    qq = qq + nr;
    pp.t = pp.t + nr;
end
A_t = sparse(A_t);
b_t = sparse((nb-1)*nr,1);
%
con_ax = [con_ax; A_t];
con_bx = [con_bx; b_t];
ctype  = [ctype, repmat('U',1,(nb-1)*nr)];
clear pp qq A_t b_t
% <---
% Linear inequalities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear equalities
% --->
%%% equalities for (b) : regularization
pp = mk_pointer_1(num);
%
A_b = zeros(nb,num.all_var);
A_b(:,pp.b+1:pp.b+nb) = eye(nb);
A_b = sparse(A_b);
b_b = ones(nb,1);
%
con_ax = [con_ax; A_b];
con_bx = [con_bx; b_b];
ctype  = [ctype, repmat('S',1,nb)];
clear pp A_b b_b

%%% equalities for (t_11) : fixed variable
pp = mk_pointer_1(num);
%
A_t_11 = sparse(1,num.all_var);
A_t_11(:,pp.t+1) = 1;
b_t_11 = 1;
%
con_ax = [con_ax; A_t_11];
con_bx = [con_bx; b_t_11];
ctype  = [ctype, repmat('S',1,1)];
clear pp A_t_11 b_t_11
% <---
% Linear equalities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quadratic inequalities
% --->
con_Qx_qc = [];
con_rx_qc = [];
con_bx_qc = [];
%%% QC inequalities: first (nb-1) inequalities
pp = mk_pointer_1(num);
%
Qqc_vabct1 = cell((nb-1)*nr,1);
rqc_vabct1 = cell((nb-1)*nr,1);
qq = 0;
for i=1:(nb-1)
    for l=1:nr
        Qqc_vabct1{qq+l} = zeros(num.all_var,num.all_var);
        Qqc_vabct1{qq+l}(pp.a+1,pp.a+1) =  data_eps(l)^2;
        Qqc_vabct1{qq+l}(pp.b+1,pp.b+1) =  data_sig(l)^2;
        Qqc_vabct1{qq+l}(pp.c+1,pp.c+1) =  1;
        Qqc_vabct1{qq+l}(pp.a+1,pp.b+1) =  data_eps(l) * data_sig(l);
        Qqc_vabct1{qq+l}(pp.b+1,pp.a+1) =  data_eps(l) * data_sig(l);
        Qqc_vabct1{qq+l}(pp.b+1,pp.c+1) = -data_sig(l);
        Qqc_vabct1{qq+l}(pp.c+1,pp.b+1) = -data_sig(l);
        Qqc_vabct1{qq+l}(pp.c+1,pp.a+1) = -data_eps(l);
        Qqc_vabct1{qq+l}(pp.a+1,pp.c+1) = -data_eps(l);
        Qqc_vabct1{qq+l} = sparse(Qqc_vabct1{qq+l});
        %
        rqc_vabct1{qq+l} = zeros(num.all_var,1);
        rqc_vabct1{qq+l}(pp.v+l) = -1;
        rqc_vabct1{qq+l}(pp.t+l)    =  param.big_M;
        rqc_vabct1{qq+l}(pp.t+nr+l) = -param.big_M;
        rqc_vabct1{qq+l} = sparse(rqc_vabct1{qq+l});
    end
    %
    qq = qq + nr;
    pp.a = pp.a + 1;
    pp.b = pp.b + 1;
    pp.c = pp.c + 1;
    pp.v = pp.v + nr;
    pp.t = pp.t + nr;
end
bqc_vabct1 = param.big_M * ones((nb-1)*nr,1);

con_Qx_qc = [con_Qx_qc; Qqc_vabct1];
con_rx_qc = [con_rx_qc; rqc_vabct1];
con_bx_qc = [con_bx_qc; bqc_vabct1];
% clear pp Qqc_vabct1 rqc_vabct1 bqc_vabct1

%%% QC inequalities: last inequality
pp = mk_pointer_1(num);
pp.v = pp.v + ((nb-1) * nr);
pp.t = pp.t + ((nb-1) * nr);
%
Qqc_vabct2 = cell(nr,1);
rqc_vabct2 = cell(nr,1);
for l=1:nr
    Qqc_vabct2{l} = zeros(num.all_var,num.all_var);
    Qqc_vabct2{l}(pp.a+nb,pp.a+nb) =  data_eps(l)^2;
    Qqc_vabct2{l}(pp.b+nb,pp.b+nb) =  data_sig(l)^2;
    Qqc_vabct2{l}(pp.c+nb,pp.c+nb) =  1;
    Qqc_vabct2{l}(pp.a+nb,pp.b+nb) =  data_eps(l) * data_sig(l);
    Qqc_vabct2{l}(pp.b+nb,pp.a+nb) =  data_eps(l) * data_sig(l);
    Qqc_vabct2{l}(pp.b+nb,pp.c+nb) = -data_sig(l);
    Qqc_vabct2{l}(pp.c+nb,pp.b+nb) = -data_sig(l);
    Qqc_vabct2{l}(pp.c+nb,pp.a+nb) = -data_eps(l);
    Qqc_vabct2{l}(pp.a+nb,pp.c+nb) = -data_eps(l);
    Qqc_vabct2{l} = sparse(Qqc_vabct2{l});
    %   
    rqc_vabct2{l} = zeros(num.all_var,1);
    rqc_vabct2{l}(pp.v+l) = -1;
    rqc_vabct2{l}(pp.t+l) =  param.big_M;
    rqc_vabct2{l} = sparse(rqc_vabct2{l});
end
bqc_vabct2 = param.big_M * ones(nr,1);

con_Qx_qc = [con_Qx_qc; Qqc_vabct2];
con_rx_qc = [con_rx_qc; rqc_vabct2];
con_bx_qc = [con_bx_qc; bqc_vabct2];
clear pp Qqc_vabct2 rqc_vabct2 bqc_vabct2
% <---
% Quadratic inequalities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make LP-file
% --->
if Flag.exe >= 2
    fprintf('       ------  create lp-file   ------ \n');
    [~] = make_qp_lp_file(filename.lp, vec_obj, con_ax, con_bx,...
        vec_lb, vec_ub, [], [], vartype, ctype,...
        con_Qx_qc, con_rx_qc, con_bx_qc);
    fprintf('       ------ lp-file completed ------ \n');
end
% <---
% Make LP-file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve IP
% --->
%%%% Solver options -->
if Flag.exe >= 2
    fid = fopen('cplex_opt.bat', 'w');
    fprintf(fid, 'set default\n');
    fprintf(fid, 'set logfile %s\n', filename.cplex_log);
    % fprintf(fid, 'set mip interval 100000\n');
    fprintf(fid, 'set mip tolerances integrality 1e-08 \n');
    fprintf(fid, 'set mip tolerances mipgap 1e-08 \n');
    fprintf(fid, 'set mip strategy miqcpstrat %g \n', Flag.miqcqstrat);
    % fprintf(fid, 'set timelimit %g\n', Flag.timelimit);
    fprintf(fid, 'read %s\n', filename.lp);
    fprintf(fid, 'optimize\n');
    fprintf(fid, 'write %s\n', filename.cplex_sol);
    fclose(fid);
end
%%%% <-- Solver options
%%%% Execute solver -->
if Flag.exe == 3
    tic;
    delete('cplex.log');
    delete(filename.cplex_log);
    delete(filename.cplex_sol);
    if Flag.display == 0
        [dos_disc.status, dos_disc.output] = dos('cplex.exe < cplex_opt.bat');
    else
        dos('cplex.exe < cplex_opt.bat');
    end
    time_mip = toc;
    fprintf('   CPU time for solving MIP = %3.3f s \n',time_mip);
end
%%%% <-- Execute solver
%%%% Read sol-file -->
if (Flag.exe == 1) || (Flag.exe == 3)
    Sol = readstruct(filename.cplex_sol,'FileType','xml');
    x_idx = zeros(num.all_var,1);
    for i=1:num.all_var
        x_name =...
            convertStringsToChars(Sol.variables.variable(i).nameAttribute);
        x_idx(i) = str2double(x_name(3:end));
    end
    x_min = zeros(num.all_var,1);
    for i=1:num.all_var
        x_min(x_idx(i)) = Sol.variables.variable(i).valueAttribute;
    end
    clear x_idx
    fprintf('   ----- read %s ----- \n',...
        filename.cplex_sol);
    [opt_v, opt_a, opt_b, opt_c, opt_t] =...
        mk_out_1(x_min,num);
    %
    fprintf('=============================================\n');
    fprintf('   #data points = %g\n', num.data );
    if sum(opt_t{nb}) == 0
        opt_nb = nb - 1;
        if sum(opt_t{nb-1}) == 0
            opt_nb = nb - 2;
            if sum(opt_t{nb-2}) == 0
                opt_nb = nb - 3;
            end
        end
    else
        opt_nb = nb;
    end
    idx_bp = zeros(opt_nb-1,1);
    for i=1:(opt_nb-1)
        idx_bp(i) = find(opt_t{i+1}, 1) - 1;
    end
    for i=1:(opt_nb-1)
        fprintf('   break point: bw/ %g & %g\n',...
        idx_bp(i), idx_bp(i)+1 );
    end
    for i=1:nb
        fprintf('   [a,b,c] = [%1.4f,%1.4f,%1.4f]\n',...
            opt_a(i), opt_b(i), opt_c(i) );
    end
    fprintf('=============================================\n');
    delete('clone*.log')
end
%%%% <-- Read sol-file
% <---
% Solve IP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalization
% --->
if (Flag.exe == 1) || (Flag.exe == 3)
    fprintf('   Normalization:\n');
    for i=1:nb
        norm_ab = norm( [opt_a(i),opt_b(i)] );
        opt_a(i) = opt_a(i) / norm_ab;
        opt_b(i) = opt_b(i) / norm_ab;
        opt_c(i) = opt_c(i) / norm_ab;
    end
    for i=1:opt_nb
        fprintf('   [a,b,c] = [%1.4f,%1.4f,%1.4f]\n',...
            opt_a(i), opt_b(i), opt_c(i) );
    end
    fprintf('=============================================\n');
end
% <---
% Normalization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare figures: end-points of line
% --->
if (Flag.exe == 1) || (Flag.exe == 3)
    intersec_pt = cell(1,opt_nb-1);
    for i=1:(opt_nb-1)
        intersec_pt{i} = zeros(2,1);
        al = opt_a(i);
        aa = opt_a(i+1);
        bt = opt_b(i);
        bb = opt_b(i+1);
        gm = opt_c(i);
        cc = opt_c(i+1);
        intersec_pt{i}(1) = ((bb*gm) - (bt*cc)) / ((al*bb) - (aa*bt));
        intersec_pt{i}(2) = ((cc*al) - (gm*aa)) / ((al*bb) - (aa*bt));
        clear al aa bt bb gm cc
    end

    delta_eps = data_eps(nr) - data_eps(1);
    x_plot = [data_eps(1) - (delta_eps/nr)];
    for i=1:(opt_nb-1)
        x_plot = [x_plot, intersec_pt{i}(1)];
    end
    x_plot = [x_plot,...
        data_eps(nr) + (delta_eps/nr)];

    y_plot =...
        [-( opt_a(1) * x_plot(1) /opt_b(1) )...
            + ( opt_c(1) / opt_b(1) )];
    for i=1:(opt_nb-1)
        y_plot = [y_plot, intersec_pt{i}(2)];
    end
    y_plot = [y_plot,...
        -( opt_a(opt_nb) * x_plot(end) /opt_b(opt_nb) )...
            + ( opt_c(opt_nb) / opt_b(opt_nb) )];

    brk_pt = zeros(opt_nb-1,2);
    for i=1:(opt_nb-1)
        brk_pt(i,1) =...
            ( data_eps(idx_bp(i)) + data_eps(idx_bp(i)+1) ) / 2;
        brk_pt(i,2) =...
            ( data_sig(idx_bp(i)) + data_sig(idx_bp(i)+1) ) / 2;
    end
end
% <---
% Prepare figures: end-points of line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures
% --->
if (Flag.exe == 1) || (Flag.exe == 3)
    figure;
    plot(data_eps, data_sig, 'go',...
        'MarkerFaceColor','w', 'MarkerSize',3, 'LineWidth',0.75);
    hold on;
    grid on;
    plot(x_plot, y_plot, 'k-', 'LineWidth',1.0);
    xlabel('Strain ($10^{-3}$ m/m)', 'Interpreter', 'latex');
    ylabel('Stress ($10^{6}$ Pa)', 'Interpreter', 'latex');
    xlim( [-6,6] );
    ylim( [-4.5,4.5] );
    axis equal;
    set(gcf,'renderer','painters');
    set(gca,'FontName','Times New Roman');
    set(gca,'FontSize',16);
end
% <---
% Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save('segmented_regress.mat',...
    'opt_a', 'opt_b', 'opt_c',...
    'opt_nb', 'intersec_pt',...
    'data_eps', 'data_sig');
