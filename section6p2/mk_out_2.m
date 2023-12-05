function [vec_ep, vec_sg, vec_u, vec_s, vec_t] = mk_out_2(x, num)

nm = num.member;
nd = num.degree;
%%%% variables : [ep]    \in (3*nm)                                %%%%
%%%% variables : [sg]    \in (3*nm)                                %%%%
%%%% variables : [u]     \in (nd)                                %%%%
%%%% variables : [s]     \in (3*nm)                                %%%%
%%%% variables : [t]     \in (3*nm)                                %%%%

nmm = 3 * nm;

vec_ep = x(1:nmm);
x = x((nmm+1):end);

vec_sg = x(1:nmm);
x = x((nmm+1):end);

vec_u = x(1:nd);
x = round( x((nd+1):end) );

vec_s = x(1:nmm);
x = round( x((nmm+1):end) );

vec_t = x(1:nmm);
x = round( x((nmm+1):end) );

if ~isempty(x)
    fprintf('\n !!!!! Error in mk_out_2.m !!!!! \n');
end

