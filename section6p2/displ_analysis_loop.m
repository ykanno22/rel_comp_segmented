clear;  
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global vec_p matL matN ref_Young ref_back intersec_pt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some parameters
% --->
param.eps  = 0.10;
param.delta = 0.10;

num_Iter = 20;

param.big_M = 10^(2);
param.bi_sec_terminate = 10^(-7);
param.d_max = 10;
param.load_factor = 24.0;

param.ca_unit = 10.0; %%%% in cm^2
param.radius_ratio = 0.90;
param.ci_coeff = 10^(-4) * (1/(4*pi))...
    * (1 + (param.radius_ratio^2))...
    / (1 - (param.radius_ratio^2));

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

data_Young  = data_sig ./ data_eps;
mean_Young  = mean(data_Young);

ref_Young = zeros(ns,1);
ref_back  = zeros(ns,1);
for i=1:ns
    ref_Young(i)  = -opt_a(i) / opt_b(i);
    ref_back(i)   =  opt_c(i) / opt_b(i);
end

filename.std = 'milp_displ_analysis';
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
% Frame data
% --->
[dll,matBt,coord_x,ir,irr,ird] = member(0);
%
nk = size(coord_x,1);   num.node   = nk;
nd = size(matBt{1},1);  num.degree = nd;
nm = size(matBt{1},2);  num.member = nm;
%
Idx_load = 8;
%
vec_cs = 10.0 * ones(nm,1);
% <---
% Frame data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load vector
% --->
vec_p0 = zeros(nd,1);
vec_p0(Idx_load)   = -1.0;
%
vec_p0 = param.load_factor * vec_p0;
% <---
% Load vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stiffness matrix
% --->
cs.area    = param.ca_unit * ones(nm,1);
cs.inertia = (param.ci_coeff * param.ca_unit) * ones(nm,1);

vec_k{1} =    mean_Young * (cs.area ./ dll);
vec_k{2} = 3* mean_Young * (cs.inertia ./ dll);
vec_k{3} =    mean_Young * (cs.inertia ./ dll);

matK =...
    (matBt{1} * sparse(diag(vec_k{1})) * matBt{1}') + ...
    (matBt{2} * sparse(diag(vec_k{2})) * matBt{2}') + ...
    (matBt{3} * sparse(diag(vec_k{3})) * matBt{3}');

ref_u = matK \ vec_p0;

matL = [
    diag(1 ./ dll) * matBt{1}';...
    diag(1 ./ dll) * matBt{2}';...
    diag(1 ./ dll) * matBt{3}'];
matN = [
    matBt{1} * diag(cs.area),...
    matBt{2} * (3 * diag(cs.inertia)),...
    matBt{3} * diag(cs.inertia)];
ref_eps = matL * ref_u;
ref_sig = mean_Young * ref_eps;
% <---
% Stiffness matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discrete optimization
% --->
fprintf('\n  ##### MILP >>>>> \n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Define variable types -->
%%%% variables : [ep]    \in (3*nm)                              %%%%
%%%% variables : [sg]    \in (3*nm)                              %%%%
%%%% variables : [u]     \in (nd)                                %%%%
%%%% variables : [s]     \in (3*nm)                              %%%%
%%%% variables : [t]     \in (3*nm)                              %%%%
nmm = 3 * nm;

char.ep = repmat('C',nmm,1);
char.sg = repmat('C',nmm,1);
char.u  = repmat('C',nd,1);
char.s  = repmat('B',nmm,1);
char.t  = repmat('B',nmm,1);

vartype =...
    [char.ep; char.sg; char.u; char.s; char.t];
clear char

num.all_var = length(vartype);
%%%% <-- Define variable types
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective function
% --->
obj.ep = repmat(0,nmm,1);
obj.sg = repmat(0,nmm,1);
obj.u  = repmat(0,nd,1);
obj.u(Idx_load) = -1;
obj.s  = repmat(0,nmm,1);
obj.t  = repmat(0,nmm,1);

vec_obj =...
    sparse([obj.ep; obj.sg; obj.u; obj.s; obj.t]);
clear obj
% <---
% Objective function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define lower/upper bounds
% --->
%%%% lower bounds
lb.ep = -Inf(nmm,1);
lb.sg = -Inf(nmm,1);
lb.u  = -Inf(nd,1);
lb.s  =  zeros(nmm,1);
lb.t  =  zeros(nmm,1);
%
vec_lb = [lb.ep; lb.sg; lb.u; lb.s; lb.t];
clear lb

%%%% upper bounds
ub.ep = Inf(nmm,1);
ub.sg = Inf(nmm,1);
ub.u  = Inf(nd,1);
ub.s  = ones(nmm,1);
ub.t  = ones(nmm,1);
%
vec_ub = [ub.ep; ub.sg; ub.u; ub.s; ub.t];
clear ub
% <---
% Define lower/upper bounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
his_load_fact = zeros(1, num_Iter+1);
his_model_u = zeros(nd, num_Iter+1);
his_min_u  = zeros(nd, num_Iter+1);
his_max_u  = zeros(nd, num_Iter+1);

incr_load_factor = 1 / num_Iter;

for Iter = 0:num_Iter
    load_fact = Iter * incr_load_factor;
    vec_p = load_fact * vec_p0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    filename.lp = strcat(filename.std,'_max.lp');
    filename.cplex_log = strcat(filename.std,'_max_cplex.log');
    filename.cplex_sol = strcat(filename.std,'_max_cplex.sol');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Linear inequalities
    % --->
    ctype = [];
    con_ax = [];
    con_bx = [];

    %%% inequalities for (s,t) : order of 0-1 variables
    pp = mk_pointer_2(num);
    %
    A_st = zeros(nmm,num.all_var);
    %
    A_st(:,(pp.s+1):(pp.s+nmm)) = -eye(nmm);
    A_st(:,(pp.t+1):(pp.t+nmm)) =  eye(nmm);
    A_st = sparse(A_st);
    b_st = sparse(nmm,1);
    %
    con_ax = [con_ax; A_st];
    con_bx = [con_bx; b_st];
    ctype  = [ctype, repmat('U',1,nmm)];
    clear pp A_st b_st

    %%% inequalities for (ep,sg,s,t) : assortment 1 (1/2)
    pp = mk_pointer_2(num);
    %
    A_epsgs = zeros(nmm,num.all_var);
    %
    A_epsgs(:,(pp.ep+1):(pp.ep+nmm)) =  cst_p(1) * eye(nmm);
    A_epsgs(:,(pp.sg+1):(pp.sg+nmm)) =  cst_q(1) * eye(nmm);
    A_epsgs(:,(pp.s+1):(pp.s+nmm))   = -param.big_M * eye(nmm);
    A_epsgs(:,(pp.t+1):(pp.t+nmm))   = -param.big_M * eye(nmm);
    A_epsgs = sparse(A_epsgs);
    b_epsgs = cst_r(1) * ones(nmm,1);
    %
    con_ax = [con_ax; A_epsgs];
    con_bx = [con_bx; b_epsgs];
    ctype  = [ctype, repmat('U',1,nmm)];
    clear pp A_epsgs b_epsgs

    %%% inequalities for (ep,sg,s,t) : assortment 1 (2/2)
    pp = mk_pointer_2(num);
    %
    A_epsgs = zeros(nmm,num.all_var);
    %
    A_epsgs(:,(pp.ep+1):(pp.ep+nmm)) = -cst_p(1) * eye(nmm);
    A_epsgs(:,(pp.sg+1):(pp.sg+nmm)) = -cst_q(1) * eye(nmm);
    A_epsgs(:,(pp.s+1):(pp.s+nmm))   =  param.big_M * eye(nmm);
    A_epsgs(:,(pp.t+1):(pp.t+nmm))   = -param.big_M * eye(nmm);
    A_epsgs = sparse(A_epsgs);
    b_epsgs = (param.big_M - cst_r(1)) * ones(nmm,1);
    %
    con_ax = [con_ax; A_epsgs];
    con_bx = [con_bx; b_epsgs];
    ctype  = [ctype, repmat('U',1,nmm)];
    clear pp A_epsgs b_epsgs

    %%% inequalities for (ep,sg,s,t) : assortment 2 (1/2)
    pp = mk_pointer_2(num);
    %
    A_epsgt = zeros(nmm,num.all_var);
    %
    A_epsgt(:,(pp.ep+1):(pp.ep+nmm)) =  cst_p(2) * eye(nmm);
    A_epsgt(:,(pp.sg+1):(pp.sg+nmm)) =  cst_q(2) * eye(nmm);
    A_epsgt(:,(pp.s+1):(pp.s+nmm))   =  param.big_M * eye(nmm);
    A_epsgt(:,(pp.t+1):(pp.t+nmm))   = -param.big_M * eye(nmm);
    A_epsgt = sparse(A_epsgt);
    b_epsgt = (param.big_M + cst_r(2)) * ones(nmm,1);
    %
    con_ax = [con_ax; A_epsgt];
    con_bx = [con_bx; b_epsgt];
    ctype  = [ctype, repmat('U',1,nmm)];
    clear pp A_epsgt b_epsgt

    %%% inequalities for (ep,sg,s,t) : assortment 2 (2/2)
    pp = mk_pointer_2(num);
    %
    A_epsgt = zeros(nmm,num.all_var);
    %
    A_epsgt(:,(pp.ep+1):(pp.ep+nmm)) = -cst_p(2) * eye(nmm);
    A_epsgt(:,(pp.sg+1):(pp.sg+nmm)) = -cst_q(2) * eye(nmm);
    A_epsgt(:,(pp.s+1):(pp.s+nmm))   =  param.big_M * eye(nmm);
    A_epsgt(:,(pp.t+1):(pp.t+nmm))   =  param.big_M * eye(nmm);
    A_epsgt = sparse(A_epsgt);
    b_epsgt = ((2*param.big_M) - cst_r(2)) * ones(nmm,1);
    %
    con_ax = [con_ax; A_epsgt];
    con_bx = [con_bx; b_epsgt];
    ctype  = [ctype, repmat('U',1,nmm)];
    clear pp A_epsgt b_epsgt

    %%% inequalities for (ep,sg,s,t) : 1st segment (1/2)
    pp = mk_pointer_2(num);
    %
    A_all = zeros(nmm,num.all_var);
    %
    A_all(:,(pp.ep+1):(pp.ep+nmm)) =  opt_a(1) * eye(nmm);
    A_all(:,(pp.sg+1):(pp.sg+nmm)) =  opt_b(1) * eye(nmm);
    A_all(:,(pp.s+1):(pp.s+nmm))   = -param.big_M * eye(nmm);
    A_all(:,(pp.t+1):(pp.t+nmm))   = -param.big_M * eye(nmm);
    A_all = sparse(A_all);
    b_all = (opt_c(1) + val_d) * ones(nmm,1);
    %
    con_ax = [con_ax; A_all];
    con_bx = [con_bx; b_all];
    ctype  = [ctype, repmat('U',1,nmm)];
    clear pp A_all b_all

    %%% inequalities for (ep,sg,s,t) : 1st segment (2/2)
    pp = mk_pointer_2(num);
    %
    A_all = zeros(nmm,num.all_var);
    %
    A_all(:,(pp.ep+1):(pp.ep+nmm)) = -opt_a(1) * eye(nmm);
    A_all(:,(pp.sg+1):(pp.sg+nmm)) = -opt_b(1) * eye(nmm);
    A_all(:,(pp.s+1):(pp.s+nmm))   = -param.big_M * eye(nmm);
    A_all(:,(pp.t+1):(pp.t+nmm))   = -param.big_M * eye(nmm);
    A_all = sparse(A_all);
    b_all = (-opt_c(1) + val_d) * ones(nmm,1);
    %
    con_ax = [con_ax; A_all];
    con_bx = [con_bx; b_all];
    ctype  = [ctype, repmat('U',1,nmm)];
    clear pp A_all b_all

    %%% inequalities for (ep,sg,s,t) : 2nd segment (1/2)
    pp = mk_pointer_2(num);
    %
    A_all = zeros(nmm,num.all_var);
    %
    A_all(:,(pp.ep+1):(pp.ep+nmm)) =  opt_a(2) * eye(nmm);
    A_all(:,(pp.sg+1):(pp.sg+nmm)) =  opt_b(2) * eye(nmm);
    A_all(:,(pp.s+1):(pp.s+nmm))   =  param.big_M * eye(nmm);
    A_all(:,(pp.t+1):(pp.t+nmm))   = -param.big_M * eye(nmm);
    A_all = sparse(A_all);
    b_all = (opt_c(2) + val_d + param.big_M) * ones(nmm,1);
    %
    con_ax = [con_ax; A_all];
    con_bx = [con_bx; b_all];
    ctype  = [ctype, repmat('U',1,nmm)];
    clear pp A_all b_all

    %%% inequalities for (ep,sg,s,t) : 2nd segment (2/2)
    pp = mk_pointer_2(num);
    %
    A_all = zeros(nmm,num.all_var);
    %
    A_all(:,(pp.ep+1):(pp.ep+nmm)) = -opt_a(2) * eye(nmm);
    A_all(:,(pp.sg+1):(pp.sg+nmm)) = -opt_b(2) * eye(nmm);
    A_all(:,(pp.s+1):(pp.s+nmm))   =  param.big_M * eye(nmm);
    A_all(:,(pp.t+1):(pp.t+nmm))   = -param.big_M * eye(nmm);
    A_all = sparse(A_all);
    b_all = (-opt_c(2) + val_d + param.big_M) * ones(nmm,1);
    %
    con_ax = [con_ax; A_all];
    con_bx = [con_bx; b_all];
    ctype  = [ctype, repmat('U',1,nmm)];
    clear pp A_all b_all

    %%% inequalities for (ep,sg,s,t) : 3rd segment (1/2)
    pp = mk_pointer_2(num);
    %
    A_all = zeros(nmm,num.all_var);
    %
    A_all(:,(pp.ep+1):(pp.ep+nmm)) =  opt_a(3) * eye(nmm);
    A_all(:,(pp.sg+1):(pp.sg+nmm)) =  opt_b(3) * eye(nmm);
    A_all(:,(pp.s+1):(pp.s+nmm))   =  param.big_M * eye(nmm);
    A_all(:,(pp.t+1):(pp.t+nmm))   =  param.big_M * eye(nmm);
    A_all = sparse(A_all);
    b_all = (opt_c(3) + val_d + (2*param.big_M)) * ones(nmm,1);
    %
    con_ax = [con_ax; A_all];
    con_bx = [con_bx; b_all];
    ctype  = [ctype, repmat('U',1,nmm)];
    clear pp A_all b_all

    %%% inequalities for (ep,sg,s,t) : 3rd segment (2/2)
    pp = mk_pointer_2(num);
    %
    A_all = zeros(nmm,num.all_var);
    %
    A_all(:,(pp.ep+1):(pp.ep+nmm)) = -opt_a(3) * eye(nmm);
    A_all(:,(pp.sg+1):(pp.sg+nmm)) = -opt_b(3) * eye(nmm);
    A_all(:,(pp.s+1):(pp.s+nmm))   =  param.big_M * eye(nmm);
    A_all(:,(pp.t+1):(pp.t+nmm))   =  param.big_M * eye(nmm);
    A_all = sparse(A_all);
    b_all = (-opt_c(3) + val_d + (2*param.big_M)) * ones(nmm,1);
    %
    con_ax = [con_ax; A_all];
    con_bx = [con_bx; b_all];
    ctype  = [ctype, repmat('U',1,nmm)];
    clear pp A_all b_all
    % <---
    % Linear inequalities
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Linear equalities
    % --->
    %%% equalities for (ep,u) : compatibility
    pp = mk_pointer_2(num);
    %
    A_epu = zeros(nmm,num.all_var);
    A_epu(:,(pp.ep+1):(pp.ep+nmm)) =  eye(nmm);
    A_epu(:,(pp.u+1):(pp.u+nd))   = -matL;
    A_epu = sparse(A_epu);
    b_epu = sparse(nmm,1);
    %
    con_ax = [con_ax; A_epu];
    con_bx = [con_bx; b_epu];
    ctype  = [ctype, repmat('S',1,nmm)];
    clear pp A_epu b_epu

    %%% equalities for (sg) : force-balance
    pp = mk_pointer_2(num);
    %
    A_sg = zeros(nd,num.all_var);
    A_sg(:,(pp.sg+1):(pp.sg+nmm)) = matN;
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
    % Make LP-file
    % --->
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
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    filename.lp = strcat(filename.std,'_min.lp');
    filename.cplex_log = strcat(filename.std,'_min_cplex.log');
    filename.cplex_sol = strcat(filename.std,'_min_cplex.sol');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Make LP-file
    % --->
    if Flag.exe >= 2
        fprintf('       ------  create lp-file   ------ \n');
        [~] = make_lp_file(filename.lp, -vec_obj, con_ax, con_bx,...
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

        fprintf('\n\n');
        fprintf('   upper bound: displ = %1.4f \n',...
            max_u(Idx_load) );
        fprintf('   lower bound: displ = %1.4f \n',...
            min_u(Idx_load) );
        fprintf('=============================================\n');
        delete('clone*.log')
    end
    %%%% <-- Read sol-file
    % <---
    % Solve IP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reference solution
    % -->
    fprintf('   Reference solution via regression data -->');
    options = optimoptions('fsolve', 'MaxIter',2000, 'MaxFunEvals',10^8,...
        'FunctionTolerance',10^(-8), 'OptimalityTolerance',10^(-8),...
        'Display','off');
    [ref_u,~,exit_fsolve] = fsolve(@ref_eq_residual, ref_u, options);
    fprintf('\n   Exit_flag of fsolve = %g\n', exit_fsolve);
    % <--
    % Reference solution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    his_ref_u(:,Iter+1)  = ref_u;
    his_min_u(:,Iter+1)  = min_u;
    his_max_u(:,Iter+1)  = max_u;
    his_load_fact(1,Iter+1) = load_fact;
end
% <---
% Discrete optimization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures
% --->
figure;
plot(his_max_u(Idx_load,:), his_load_fact, 'r<-',...
    'LineWidth', 1.0, 'MarkerSize',6);
hold on;
plot(his_min_u(Idx_load,:), his_load_fact, 'r>-',...
    'LineWidth', 1.0, 'MarkerSize',6);
plot(his_ref_u(Idx_load,:), his_load_fact, 'bs:',...
    'LineWidth', 1.0, 'MarkerSize',6);
xlabel('Displacement (mm)', 'Interpreter', 'latex');
ylabel('Load factor, $\lambda$', 'Interpreter', 'latex');
set(gcf,'renderer','painters');
set(gca,'FontName','Times New Roman');
set(gca,'FontSize',16);
% <---
% Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

