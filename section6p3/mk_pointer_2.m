function [pp] = mk_pointer_2(num)

nm = num.member;
nd = num.degree;
%%%% variables : [ep]    \in (nm)                                %%%%
%%%% variables : [sg]    \in (nm)                                %%%%
%%%% variables : [u]     \in (nd)                                %%%%
%%%% variables : [s]     \in (nm)                                %%%%
%%%% variables : [t]     \in (nm)                                %%%%
%
pp.ep = 0;
pp.sg = pp.ep + nm;
pp.u  = pp.sg + nm;
pp.s  = pp.u  + nd;
pp.t  = pp.s  + nm;
