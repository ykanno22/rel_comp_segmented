function [pp] = mk_pointer_2(num)

nm = num.member;
nd = num.degree;
%%%% variables : [ep]    \in (3*nm)                              %%%%
%%%% variables : [sg]    \in (3*nm)                              %%%%
%%%% variables : [u]     \in (nd)                                %%%%
%%%% variables : [s]     \in (3*nm)                              %%%%
%%%% variables : [t]     \in (3*nm)                              %%%%
%
pp.ep = 0;
pp.sg = pp.ep + (3*nm);
pp.u  = pp.sg + (3*nm);
pp.s  = pp.u  + nd;
pp.t  = pp.s  + (3*nm);
