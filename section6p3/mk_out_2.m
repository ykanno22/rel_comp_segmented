function [vec_ep, vec_sg, vec_u, vec_s, vec_t] = mk_out_2(x, num)

nm = num.member;
nd = num.degree;
%%%% variables : [ep]    \in (nm)                                %%%%
%%%% variables : [sg]    \in (nm)                                %%%%
%%%% variables : [u]     \in (nd)                                %%%%
%%%% variables : [s]     \in (nm)                                %%%%
%%%% variables : [t]     \in (nm)                                %%%%

vec_ep = x(1:nm);
x = x((nm+1):end);

vec_sg = x(1:nm);
x = x((nm+1):end);

vec_u = x(1:nd);
x = round( x((nd+1):end) );

vec_s = x(1:nm);
x = round( x((nm+1):end) );

vec_t = x(1:nm);
x = round( x((nm+1):end) );

if ~isempty(x)
    fprintf('\n !!!!! Error in mk_out_2.m !!!!! \n');
end

