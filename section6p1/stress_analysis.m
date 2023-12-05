clear;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global vec_p matL matN ref_Young ref_back intersec_pt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some parameters
% --->
param.eps  = 0.10;
param.delta = 0.10;

param.big_M = 10^(2);
param.bi_sec_terminate = 10^(-7);
param.d_max = 10;
param.load_factor = 21.0;

% Flag.exe = 1; %%%% read sol-file and create figure
% Flag.exe = 2; %%%% create lp-file
Flag.exe = 3; %%%% solve problem
% Flag.miqcqstrat = 1;  %%% QP-base
Flag.miqcqstrat = 2;  %%% LP-base
% Flag.display = 1;
Flag.display = 0;

load('segmented_regress.mat');

nr = length(data_eps);
ns = opt_nb;  %%% #segments
clear opt_nb
 
num.data = nr;
num.segment = ns;

opt_a = opt_a(1:ns);
opt_b = opt_b(1:ns);
opt_c = opt_c(1:ns);

ref_Young = zeros(ns,1);
ref_back  = zeros(ns,1);
for i=1:ns
    ref_Young(i)  = -opt_a(i) / opt_b(i);
    ref_back(i)   =  opt_c(i) / opt_b(i);
end

filename.std = 'milp_stress_analysis';
% <---
% Some parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Borderlines
% --->
for i=1:(ns-1)
    cst_p(i) = opt_a(i) - opt_a(i+1);
    cst_q(i) = opt_b(i) - opt_b(i+1);
    cst_r(i) = (cst_p(i) * intersec_pt{i}(1))...
        + (cst_q(i) * intersec_pt{i}(2));
end
for i=1:(ns-1)
    if cst_p(i) < 0
        cst_p(i) = -cst_p(i);
        cst_q(i) = -cst_q(i);
        cst_r(i) = -cst_r(i);
    end
end
% <---
% Borderlines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assortment
% --->
SetS = cell(1,ns);
SetS{1} = [];
SetS{2} = [];
SetS{3} = [];
for l=1:num.data
    xx = data_eps(l);
    yy = data_sig(l);
    f_val(1) = (cst_p(1) * xx) + (cst_q(1) * yy) - cst_r(1);
    f_val(2) = (cst_p(2) * xx) + (cst_q(2) * yy) - cst_r(2);
    if f_val(1) <= 0
        SetS{1} = [SetS{1}, l];
    elseif f_val(2) <= 0
        SetS{2} = [SetS{2}, l];
    else
        SetS{3} = [SetS{3}, l];
    end
end
clear xx yy f_val
% <---
% Assortment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% #points included in uncertainty set
% --->
num_confi = num.data;
cur_delta = (1-param.eps)^num.data;

if cur_delta > param.delta
    fprintf(' confidence level "delta" is too small \n');
    return;
end

while cur_delta <= param.delta
    num_confi = num_confi - 1;
    cur_delta = cur_delta + ...
        ( nchoosek(num.data,num_confi)...
        * ( (1-param.eps)^num_confi )...
        * ( param.eps^(num.data-num_confi) ) );
end

fprintf(' ==================================================== \n');
fprintf('   "1 - %4.3f reliability" with "1 - %4.3f confidence" \n',...
    param.eps, param.delta );
fprintf('   (#conf-1)/(#data) = %g/%g: LHS = %4.3e >  delta = %4.3e \n',...
    num_confi, num.data, cur_delta, param.delta );
num_confi = num_confi + 1;
cur_delta = 0;
for jj=num_confi:num.data
    cur_delta = cur_delta +...
        ( nchoosek(num.data,jj)...
        * ( (1-param.eps)^jj )...
        * ( param.eps^(num.data-jj) ) );
end
fprintf('   (#conf)  /(#data) = %g/%g: LHS = %4.3e <= delta = %4.3e \n',...
    num_confi, num.data, cur_delta, param.delta );
fprintf(' ==================================================== \n');

clear cur_delta jj
% <---
% #points included in uncertainty set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bi-section method for val_d
% --->
param.d_min = 0;
val_d = param.d_max / 2;
delta_d = param.d_max - param.d_min;

while delta_d > param.bi_sec_terminate
    num_inluded = 0;
    for i=1:ns
        for l=SetS{i}
            xx = data_eps(l);
            yy = data_sig(l);
            LHS = abs(  (opt_a(i) * xx) + (opt_b(i) * yy) - opt_c(i)  );
            if LHS <= val_d
                num_inluded = num_inluded + 1;
            end
        end
    end
    if num_inluded >= num_confi
        param.d_max = val_d;
    else
        param.d_min = val_d;
    end
    delta_d = param.d_max - param.d_min;
    val_d = (param.d_max + param.d_min) / 2;
end
for i=1:ns
    fprintf('   [a,b,c] = [%1.4f,%1.4f,%1.4f]\n',...
        opt_a(i), opt_b(i), opt_c(i) );
end
fprintf('   d = %4.5e;   thr_for_bi-section: %4.3e\n',...
    val_d, param.bi_sec_terminate );
fprintf(' ==================================================== \n');
clear xx yy LHS num_inluded delta_d
% <---
% Bi-section method for val_d
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Truss data
% --->
nx = 3;
ny = 2;
[dll,matBt,coord_x,ir,irr,ird] = member(nx,ny);
%
nk = size(coord_x,1);  num.node   = nk;
nd = size(matBt,1);    num.degree = nd;
nm = size(matBt,2);    num.member = nm;
%
vec_cs = 10.0 * ones(nm,1);
% <---
% Truss data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load vector
% --->
vec_p = zeros(nd,1);
vec_p(end)   = -1.0;
vec_p(end-2) = -1.0;
%
vec_p = param.load_factor * vec_p;
% <---
% Load vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stiffness matrix & reference solution
% --->
matK = matBt * sparse(diag(ref_Young(2) * vec_cs ./ dll)) * matBt';
ref_u = matK \ vec_p;
matL  = diag(1 ./ dll) * matBt';
matN  = matBt * diag(vec_cs);
%
fprintf('   Reference solution via regression data -->');
options = optimoptions('fsolve', 'MaxIter',2000, 'MaxFunEvals',10^8,...
    'FunctionTolerance',10^(-8), 'OptimalityTolerance',10^(-8),...
    'Display','iter');
[ref_u,~,exit_fsolve] = fsolve(@ref_eq_residual, ref_u, options);
ref_eps = matL * ref_u;
ref_sig = zeros(nm,1);
for i=1:nm
    strain = ref_eps(i);
    if strain <= intersec_pt{1}(1)
        ref_sig(i) = (ref_Young(1) * strain) + ref_back(1);
    elseif strain >= intersec_pt{2}(1)
        ref_sig(i) = (ref_Young(3) * strain) + ref_back(3);
    else
        ref_sig(i) = (ref_Young(2) * strain) + ref_back(2);
    end
end
fprintf('\n   Exit_flag of fsolve = %g\n', exit_fsolve);
fprintf('   strain & stress (reference solution)\n');
for i=1:nm
    fprintf('      member %g: (%1.4f, %1.4f)\n',...
        i, ref_eps(i), ref_sig(i));
end
fprintf(' ==================================================== \n');
% <---
% Stiffness matrix & reference solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discrete optimization
% --->
fprintf('\n  ##### MILP >>>>> \n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename.lp = strcat(filename.std,'_max.lp');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Define variable types -->
%%%% variables : [ep]    \in (nm)                                %%%%
%%%% variables : [sg]    \in (nm)                                %%%%
%%%% variables : [u]     \in (nd)                                %%%%
%%%% variables : [s]     \in (nm)                                %%%%
%%%% variables : [t]     \in (nm)                                %%%%
char.ep = repmat('C',nm,1);
char.sg = repmat('C',nm,1);
char.u  = repmat('C',nd,1);
char.s  = repmat('B',nm,1);
char.t  = repmat('B',nm,1);

vartype =...
    [char.ep; char.sg; char.u; char.s; char.t];
clear char

num.all_var = length(vartype);
%%%
% <-- Define variable types
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define lower/upper bounds
% --->
%%%% lower bounds
lb.ep = -Inf(nm,1);
lb.sg = -Inf(nm,1);
lb.u  = -Inf(nd,1);
lb.s  =  zeros(nm,1);
lb.t  =  zeros(nm,1);
%
vec_lb = [lb.ep; lb.sg; lb.u; lb.s; lb.t];
clear lb

%%%% upper bounds
ub.ep = Inf(nm,1);
ub.sg = Inf(nm,1);
ub.u  = Inf(nd,1);
ub.s  = ones(nm,1);
ub.t  = ones(nm,1);
%
vec_ub = [ub.ep; ub.sg; ub.u; ub.s; ub.t];
clear ub
% <---
% Define lower/upper bounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear inequalities
% --->
ctype = [];
con_ax = [];
con_bx = [];

%%% inequalities for (s,t) : order of 0-1 variables
pp = mk_pointer_2(num);
%
A_st = zeros(nm,num.all_var);
%
A_st(:,(pp.s+1):(pp.s+nm)) = -eye(nm);
A_st(:,(pp.t+1):(pp.t+nm)) =  eye(nm);
A_st = sparse(A_st);
b_st = sparse(nm,1);
%
con_ax = [con_ax; A_st];
con_bx = [con_bx; b_st];
ctype  = [ctype, repmat('U',1,nm)];
clear pp A_st b_st

%%% inequalities for (ep,sg,s,t) : assortment 1 (1/2)
pp = mk_pointer_2(num);
%
A_epsgs = zeros(nm,num.all_var);
%
A_epsgs(:,(pp.ep+1):(pp.ep+nm)) =  cst_p(1) * eye(nm);
A_epsgs(:,(pp.sg+1):(pp.sg+nm)) =  cst_q(1) * eye(nm);
A_epsgs(:,(pp.s+1):(pp.s+nm))   = -param.big_M * eye(nm);
A_epsgs(:,(pp.t+1):(pp.t+nm))   = -param.big_M * eye(nm);
A_epsgs = sparse(A_epsgs);
b_epsgs = cst_r(1) * ones(nm,1);
%
con_ax = [con_ax; A_epsgs];
con_bx = [con_bx; b_epsgs];
ctype  = [ctype, repmat('U',1,nm)];
clear pp A_epsgs b_epsgs

%%% inequalities for (ep,sg,s,t) : assortment 1 (2/2)
pp = mk_pointer_2(num);
%
A_epsgs = zeros(nm,num.all_var);
%
A_epsgs(:,(pp.ep+1):(pp.ep+nm)) = -cst_p(1) * eye(nm);
A_epsgs(:,(pp.sg+1):(pp.sg+nm)) = -cst_q(1) * eye(nm);
A_epsgs(:,(pp.s+1):(pp.s+nm))   =  param.big_M * eye(nm);
A_epsgs(:,(pp.t+1):(pp.t+nm))   = -param.big_M * eye(nm);
A_epsgs = sparse(A_epsgs);
b_epsgs = (param.big_M - cst_r(1)) * ones(nm,1);
%
con_ax = [con_ax; A_epsgs];
con_bx = [con_bx; b_epsgs];
ctype  = [ctype, repmat('U',1,nm)];
clear pp A_epsgs b_epsgs

%%% inequalities for (ep,sg,s,t) : assortment 2 (1/2)
pp = mk_pointer_2(num);
%
A_epsgt = zeros(nm,num.all_var);
%
A_epsgt(:,(pp.ep+1):(pp.ep+nm)) =  cst_p(2) * eye(nm);
A_epsgt(:,(pp.sg+1):(pp.sg+nm)) =  cst_q(2) * eye(nm);
A_epsgt(:,(pp.s+1):(pp.s+nm))   =  param.big_M * eye(nm);
A_epsgt(:,(pp.t+1):(pp.t+nm))   = -param.big_M * eye(nm);
A_epsgt = sparse(A_epsgt);
b_epsgt = (param.big_M + cst_r(2)) * ones(nm,1);
%
con_ax = [con_ax; A_epsgt];
con_bx = [con_bx; b_epsgt];
ctype  = [ctype, repmat('U',1,nm)];
clear pp A_epsgt b_epsgt

%%% inequalities for (ep,sg,s,t) : assortment 2 (2/2)
pp = mk_pointer_2(num);
%
A_epsgt = zeros(nm,num.all_var);
%
A_epsgt(:,(pp.ep+1):(pp.ep+nm)) = -cst_p(2) * eye(nm);
A_epsgt(:,(pp.sg+1):(pp.sg+nm)) = -cst_q(2) * eye(nm);
A_epsgt(:,(pp.s+1):(pp.s+nm))   =  param.big_M * eye(nm);
A_epsgt(:,(pp.t+1):(pp.t+nm))   =  param.big_M * eye(nm);
A_epsgt = sparse(A_epsgt);
b_epsgt = ((2*param.big_M) - cst_r(2)) * ones(nm,1);
%
con_ax = [con_ax; A_epsgt];
con_bx = [con_bx; b_epsgt];
ctype  = [ctype, repmat('U',1,nm)];
clear pp A_epsgt b_epsgt

%%% inequalities for (ep,sg,s,t) : 1st segment (1/2)
pp = mk_pointer_2(num);
%
A_all = zeros(nm,num.all_var);
%
A_all(:,(pp.ep+1):(pp.ep+nm)) =  opt_a(1) * eye(nm); 
A_all(:,(pp.sg+1):(pp.sg+nm)) =  opt_b(1) * eye(nm);
A_all(:,(pp.s+1):(pp.s+nm))   = -param.big_M * eye(nm);
A_all(:,(pp.t+1):(pp.t+nm))   = -param.big_M * eye(nm);
A_all = sparse(A_all);
b_all = (opt_c(1) + val_d) * ones(nm,1);
%
con_ax = [con_ax; A_all];
con_bx = [con_bx; b_all];
ctype  = [ctype, repmat('U',1,nm)];
clear pp A_all b_all

%%% inequalities for (ep,sg,s,t) : 1st segment (2/2)
pp = mk_pointer_2(num);
%
A_all = zeros(nm,num.all_var);
%
A_all(:,(pp.ep+1):(pp.ep+nm)) = -opt_a(1) * eye(nm); 
A_all(:,(pp.sg+1):(pp.sg+nm)) = -opt_b(1) * eye(nm);
A_all(:,(pp.s+1):(pp.s+nm))   = -param.big_M * eye(nm);
A_all(:,(pp.t+1):(pp.t+nm))   = -param.big_M * eye(nm);
A_all = sparse(A_all);
b_all = (-opt_c(1) + val_d) * ones(nm,1);
%
con_ax = [con_ax; A_all];
con_bx = [con_bx; b_all];
ctype  = [ctype, repmat('U',1,nm)];
clear pp A_all b_all

%%% inequalities for (ep,sg,s,t) : 2nd segment (1/2)
pp = mk_pointer_2(num);
%
A_all = zeros(nm,num.all_var);
%
A_all(:,(pp.ep+1):(pp.ep+nm)) =  opt_a(2) * eye(nm); 
A_all(:,(pp.sg+1):(pp.sg+nm)) =  opt_b(2) * eye(nm);
A_all(:,(pp.s+1):(pp.s+nm))   =  param.big_M * eye(nm);
A_all(:,(pp.t+1):(pp.t+nm))   = -param.big_M * eye(nm);
A_all = sparse(A_all);
b_all = (opt_c(2) + val_d + param.big_M) * ones(nm,1);
%
con_ax = [con_ax; A_all];
con_bx = [con_bx; b_all];
ctype  = [ctype, repmat('U',1,nm)];
clear pp A_all b_all

%%% inequalities for (ep,sg,s,t) : 2nd segment (2/2)
pp = mk_pointer_2(num);
%
A_all = zeros(nm,num.all_var);
%
A_all(:,(pp.ep+1):(pp.ep+nm)) = -opt_a(2) * eye(nm); 
A_all(:,(pp.sg+1):(pp.sg+nm)) = -opt_b(2) * eye(nm);
A_all(:,(pp.s+1):(pp.s+nm))   =  param.big_M * eye(nm);
A_all(:,(pp.t+1):(pp.t+nm))   = -param.big_M * eye(nm);
A_all = sparse(A_all);
b_all = (-opt_c(2) + val_d + param.big_M) * ones(nm,1);
%
con_ax = [con_ax; A_all];
con_bx = [con_bx; b_all];
ctype  = [ctype, repmat('U',1,nm)];
clear pp A_all b_all

%%% inequalities for (ep,sg,s,t) : 3rd segment (1/2)
pp = mk_pointer_2(num);
%
A_all = zeros(nm,num.all_var);
%
A_all(:,(pp.ep+1):(pp.ep+nm)) =  opt_a(3) * eye(nm); 
A_all(:,(pp.sg+1):(pp.sg+nm)) =  opt_b(3) * eye(nm);
A_all(:,(pp.s+1):(pp.s+nm))   =  param.big_M * eye(nm);
A_all(:,(pp.t+1):(pp.t+nm))   =  param.big_M * eye(nm);
A_all = sparse(A_all);
b_all = (opt_c(3) + val_d + (2*param.big_M)) * ones(nm,1);
%
con_ax = [con_ax; A_all];
con_bx = [con_bx; b_all];
ctype  = [ctype, repmat('U',1,nm)];
clear pp A_all b_all

%%% inequalities for (ep,sg,s,t) : 3rd segment (2/2)
pp = mk_pointer_2(num);
%
A_all = zeros(nm,num.all_var);
%
A_all(:,(pp.ep+1):(pp.ep+nm)) = -opt_a(3) * eye(nm); 
A_all(:,(pp.sg+1):(pp.sg+nm)) = -opt_b(3) * eye(nm);
A_all(:,(pp.s+1):(pp.s+nm))   =  param.big_M * eye(nm);
A_all(:,(pp.t+1):(pp.t+nm))   =  param.big_M * eye(nm);
A_all = sparse(A_all);
b_all = (-opt_c(3) + val_d + (2*param.big_M)) * ones(nm,1);
%
con_ax = [con_ax; A_all];
con_bx = [con_bx; b_all];
ctype  = [ctype, repmat('U',1,nm)];
clear pp A_all b_all
% <---
% Linear inequalities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear equalities
% --->
%%% equalities for (ep,u) : compatibility
pp = mk_pointer_2(num);
%
A_epu = zeros(nm,num.all_var);
A_epu(:,(pp.ep+1):(pp.ep+nm)) =  eye(nm);
A_epu(:,(pp.u+1):(pp.u+nd))   = -matL;
A_epu = sparse(A_epu);
b_epu = sparse(nm,1);
%
con_ax = [con_ax; A_epu];
con_bx = [con_bx; b_epu];
ctype  = [ctype, repmat('S',1,nm)];
clear pp A_epu b_epu

%%% equalities for (sg) : force-balance
pp = mk_pointer_2(num);
%
A_sg = zeros(nd,num.all_var);
A_sg(:,(pp.sg+1):(pp.sg+nm)) = matN;
A_sg = sparse(A_sg);
b_sg = sparse(vec_p);
%
con_ax = [con_ax; A_sg];
con_bx = [con_bx; b_sg];
ctype  = [ctype, repmat('S',1,nd)];
clear pp A_sg b_sg
% <---
% Linear equalities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
max_sig_all = zeros(nm,1);
min_sig_all = zeros(nm,1);
max_eps_each = cell(nm,1);
max_sig_each = cell(nm,1);
min_eps_each = cell(nm,1);
min_sig_each = cell(nm,1);
%
for iMember = 1:nm
    filename.cplex_log =...
        strcat(filename.std,num2str(iMember),'_max_cplex.log');
    filename.cplex_sol =...
        strcat(filename.std,num2str(iMember),'_max_cplex.sol');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Make LP-file
    % --->
    pp = mk_pointer_2(num);
    vec_obj = zeros(num.all_var,1);
    vec_obj(pp.sg+iMember) = -1;
    vec_obj = sparse(vec_obj);

    if Flag.exe >= 2
        fprintf('       ------  create lp-file   ------ \n');
        [~] = make_lp_file(filename.lp, vec_obj, con_ax, con_bx,...
            vec_lb, vec_ub, vartype, ctype);
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
            [dos_IP.status, dos_IP.output] = dos('cplex.exe < cplex_opt.bat');
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
        clear x_idx Sol
        fprintf('   ----- read %s ----- \n', filename.cplex_sol);
        [max_eps, max_sig, max_u, max_s, max_t] = mk_out_2(x_min, num);

        fprintf('\n');
        fprintf('   member %g: stress upper bound = %1.4f \n',...
            iMember, max_sig(iMember) );
        fprintf('\n');
        fprintf('=============================================\n');
    end
    %%%% <-- Read sol-file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    max_sig_all(iMember) = max_sig(iMember);
    max_eps_each{iMember} = max_eps;
    max_sig_each{iMember} = max_sig;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
for iMember = 1:nm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    filename.lp = strcat(filename.std,'_min.lp');
    filename.cplex_log =...
        strcat(filename.std,num2str(iMember),'_min_cplex.log');
    filename.cplex_sol =...
        strcat(filename.std,num2str(iMember),'_min_cplex.sol');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Make LP-file
    % --->
    pp = mk_pointer_2(num);
    vec_obj = zeros(num.all_var,1);
    vec_obj(pp.sg+iMember) = 1;
    vec_obj = sparse(vec_obj);

    if Flag.exe >= 2
        fprintf('       ------  create lp-file   ------ \n');
        [~] = make_lp_file(filename.lp, vec_obj, con_ax, con_bx,...
            vec_lb, vec_ub, vartype, ctype);
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
            [dos_IP.status, dos_IP.output] = dos('cplex.exe < cplex_opt.bat');
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
        clear x_idx Sol
        fprintf('   ----- read %s ----- \n', filename.cplex_sol);
        [min_eps, min_sig, min_u, min_s, min_t] = mk_out_2(x_min, num);

        fprintf('\n');
        fprintf('   member %g: stress lower bound = %1.4f \n',...
            iMember, min_sig(iMember) );
        fprintf('\n');
        fprintf('=============================================\n');
    end
    %%%% <-- Read sol-file
    % <---
    % Solve IP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    min_sig_all(iMember) = min_sig(iMember);
    min_eps_each{iMember} = min_eps;
    min_sig_each{iMember} = min_sig;
end
fprintf('    load factor = %1.2f[kN]\n', 0.1 * param.load_factor);
fprintf('    eps = %1.4f; delta = %1.4f\n', param.eps, param.delta);
fprintf('=============================================\n');
delete('clone*.log')
delete('cplex.log');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% <---
% Discrete optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare figures: end-points of line
% --->
intersec_pts = cell(2,ns-1);
for i=1:(ns-1)
    intersec_pts{1,i} = zeros(2,1);
    al = opt_a(i);
    aa = opt_a(i+1);
    bt = opt_b(i);
    bb = opt_b(i+1);
    gm = opt_c(i);
    cc = opt_c(i+1);
    Numer = (bb*gm) - (bt*cc) + ((bt-bb) * val_d);
    intersec_pts{1,i}(1) = Numer / ((al*bb) - (aa*bt));
    Numer = (cc*al) - (gm*aa) + ((aa-al) * val_d);
    intersec_pts{1,i}(2) = Numer / ((al*bb) - (aa*bt));
    clear al aa bt bb gm cc Numer
end
for i=1:(ns-1)
    intersec_pts{2,i} = zeros(2,1);
    al = opt_a(i);
    aa = opt_a(i+1);
    bt = opt_b(i);
    bb = opt_b(i+1);
    gm = opt_c(i);
    cc = opt_c(i+1);
    Numer = (bb*gm) - (bt*cc) + ((bb-bt) * val_d);
    intersec_pts{2,i}(1) = Numer / ((al*bb) - (aa*bt));
    Numer = (cc*al) - (gm*aa) + ((al-aa) * val_d);
    intersec_pts{2,i}(2) = Numer / ((al*bb) - (aa*bt));
    clear al aa bt bb gm cc Numer
end

delta_eps = data_eps(nr) - data_eps(1);

x_plot = cell(1,ns-1);
for i=(1:ns-1)
    x_plot{i} = [data_eps(1) - (delta_eps/nr),...
        intersec_pts{i,1}(1), intersec_pts{i,2}(1),...
        data_eps(nr) + (delta_eps/nr)];
end
clear delta_eps

y_plot = cell(1,ns-1);
for i=(1:ns-1)
    y_plot{i} =...
        [-( opt_a(1) * x_plot{i}(1) /opt_b(1) )...
        + ( (opt_c(1) + ((-1)^(i)*val_d)) / opt_b(1) ),...
        intersec_pts{i,1}(2), intersec_pts{i,2}(2),...
        -( opt_a(ns) * x_plot{i}(end) /opt_b(ns) )...
        + ( (opt_c(ns) + ((-1)^(i)*val_d)) / opt_b(ns) )];
end
% <---
% Prepare figures: end-points of line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures
% --->
figure;
plot(data_eps, data_sig, 'go',...
    'MarkerFaceColor','w', 'MarkerSize',3);
hold on;
grid on;
plot(x_plot{1}, y_plot{1}, 'k-', 'LineWidth',0.5);
plot(x_plot{2}, y_plot{2}, 'k-', 'LineWidth',0.5);
plot(ref_eps, ref_sig, 'bs',...
    'MarkerSize',4, 'MarkerFaceColor','b', 'LineWidth',0.75);
if (Flag.exe == 1) || (Flag.exe == 3)
    for i=1:nm
        plot(max_eps_each{1}, max_sig_each{1}, 'r<',...
            'MarkerSize',4, 'MarkerFaceColor','w', 'LineWidth',0.75);
        plot(min_eps_each{1}, min_sig_each{1}, 'r>',...
            'MarkerSize',4, 'MarkerFaceColor','w', 'LineWidth',0.75);
    end
end
xlabel('Strain ($10^{-3}$ m/m)', 'Interpreter', 'latex');
ylabel('Stress ($10^{6}$ Pa)', 'Interpreter', 'latex');
xlim( [-6,6] );
ylim( [-4.5,4.5] );
axis equal;
set(gcf,'renderer','painters');
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',16);

figure;
hold on;
plot((1:nm), ref_sig, 'bs',...
    'MarkerSize',4, 'MarkerFaceColor','b', 'LineWidth',0.75);
plot((1:nm), max_sig_all, 'rv',...
    'MarkerSize',6, 'MarkerFaceColor','w', 'LineWidth',0.75);
plot((1:nm), min_sig_all, 'r^',...
    'MarkerSize',6, 'MarkerFaceColor','w', 'LineWidth',0.75);
xlim( [0,(nm+1)] );
ylim( [1.1*min(min_sig_all), 1.1*max(max_sig_all)] );
xlabel('Member index', 'Interpreter', 'latex');
ylabel('Stress ($10^{6}$ Pa)', 'Interpreter', 'latex');
set(gcf,'renderer','painters');
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',16);
% <---
% Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
