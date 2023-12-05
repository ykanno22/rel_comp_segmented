function [vec_v, vec_a, vec_b, vec_c, vec_t] = mk_out_1(x, num)

nr = num.data;
nb = num.b_point;
%%%% variables : [v]      \in (nb*nr)                              %%%%
%%%% variables : [a]      \in (nb)                                 %%%%
%%%% variables : [b]      \in (nb)                                 %%%%
%%%% variables : [c]      \in (nb)                                 %%%%
%%%% variables : [t]      \in (nb*nr)                              %%%%

vec_v_tmp = x(1:nb*nr);
x = x((nb*nr)+1:end);
vec_v = cell(1,nb);
for p=1:nb
    vec_v{p}  = vec_v_tmp(1:nr);
    vec_v_tmp = vec_v_tmp((nr+1):end);
end
clear vec_v_tmp

vec_a = x(1:nb);
x = x(nb+1:end);

vec_b = x(1:nb);
x = x(nb+1:end);

vec_c = x(1:nb);
x = x(nb+1:end);

vec_t_tmp = round( x(1:nb*nr) );
x = x((nb*nr)+1:end);
vec_t = cell(1,nb);
for p=1:nb
    vec_t{p}  = vec_t_tmp(1:nr);
    vec_t_tmp = vec_t_tmp((nr+1):end);
end
clear vec_v_tmp

if ~isempty(x)
    fprintf('\n !!!!! Error in mk_out_1.m !!!!! \n');
end

