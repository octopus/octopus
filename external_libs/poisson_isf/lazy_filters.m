(*Produces the interpolating-2m wavelet filter*)
m=7;
p[i_,y_]:=Product[If[j==i,1,(y-j)/(i-j)],{j,-m+1,m}];

h[q_]:=If[(Mod[q,2]==0)||(q>2m),0,p[-Quotient[q,2],1/2]];
h[0]:=1;

ht=Table[h[Abs[s]],{s,1-2m,2m-1}];
FortranForm[N[ht,18]]

