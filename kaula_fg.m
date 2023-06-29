clear
clc
close all


l=input ('l = ');
m=input ('m = ');
p=input ('p = ');
q=input ('q = ');

f=ffun(l,m,p);
g=gfun(l,p,q);


fprintf('F = %s \n',char(f))
fprintf('G = %s \n',char(g))



function f=ffun(l,m,p)

syms inc t s c

k=fix((l-m)/2);



f=symsum(((factorial(2*l-2*t))/(factorial(t)*factorial(l-t)*factorial(l-m-2*t)*2^(2*l-2*t))*((sin(inc)^(l-m-2*t))))*symsum(nchoosek(m,s)*((cos(inc))^s)*symsum(nchoosek(l-m-2*t+s,c)*nchoosek(m-s,p-t-c)*(-1)^(c-k),c,0,5),s,0,m),t,0,min(p,k));




end



function g=gfun(l,p,q)

syms e k rq rp d s t

if p>l/2
    pp=l-p;
    qp=-q;
else
    pp=p;
    qp=q;
end


if qp<0
    hp=k;
    hq=k-qp;
else
    hp=k+qp;
    hq=k;
end



beta=e/(1+sqrt(1-e^2));

g=(-1)^abs(q)*(1+beta^2)^l*beta^abs(q)*symsum(symsum(nchoosek(2*pp-2*l,hp-rp)*(((-1)^rp)/factorial(rp))*((l-2*pp+qp)*e/(2*beta))^rp,rp,0,hp)*symsum(nchoosek(-2*pp,hq-rq)/factorial(rq)*((l-2*pp+qp)*e/(2*beta))^rq,rq,0,hq)*beta^(2*k),k,0,10);


end
