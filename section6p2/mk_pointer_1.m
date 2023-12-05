function [pp] = mk_pointer_1(num)

nr = num.data;
nb = num.b_point;
%%%% variables : [v]      \in (nb*nr)                              %%%%
%%%% variables : [a]      \in (nb)                                 %%%%
%%%% variables : [b]      \in (nb)                                 %%%%
%%%% variables : [c]      \in (nb)                                 %%%%
%%%% variables : [t]      \in (nr)                                 %%%%
%
pp.v  = 0;
pp.a  = pp.v  + (nb*nr);
pp.b  = pp.a  + nb;
pp.c  = pp.b  + nb;
pp.t  = pp.c  + nb;
