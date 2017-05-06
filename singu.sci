L=10;
x=0:0.001:L;
F=(91/3)*10^3;
W=24375;
M=2.75;
S=3.8;
I=462*10^-6;
E=200*10^9;
K1=1*10^6;
K2=1*10^6;

A=[1,1,1,0;
   0,S,L,1;
   (S^3)/6,((E*I)/K1),0,-(S^2)/2;
   (L^3)/6,((L-S)^3)/6,((E*I)/K2),-(L^2)/2];

b=[-F-W*(L-S);
    F*(-1-M)-W*(L-S)*(((L-S)/2)+S);
    F*(((S-M)^2/2)-(S-M)^3/6);
    F*(((L-M)^2/2)-(L-M)^3/6)-W*(L-S)^4/24];

[x0,nsA]=linsolve(A,b);

Fa=x0(1,1);
Fc=x0(2,1);
Fd=x0(3,1);
Mr=x0(4,1);

Vx=Fa*singu0(x,0)-F*singu0(x,M)+Fc*singu0(x,S)-W*singu1(x,S)+W*singu1(x,L)+Fd*singu0(x,L);
Mx=Fa*singu1(x,0)-F*singu1(x,M)+Fc*singu1(x,S)-(W/2)*singu2(x,S)+(W/2)*singu2(x,L)+Fd*singu1(x,L)-Mr*singu0(x,0)+F*singu0(x,M);

function c=singu0(x,a)
    c=x>=a;
endfunction

function c=singu1(x,a)
    c=(x-a).*(x>=a);
endfunction

function c=singu2(x,a)
    c=((x-a).^2).*(x>=a);
endfunction

function dibujar()
    subplot(2,1,1)
    plot2d(x,Vx);
    a=gca(); // Handle on axes entity
    a.x_location = "origin";
    a.y_location = "origin";
    e = gce();
    e.children.thickness = 3;
    subplot(2,1,2)
    plot(x,Mx);
    a=gca(); // Handle on axes entity
    a.x_location = "origin";
    a.y_location = "origin";
    e = gce();
    e.children.thickness = 3;
endfunction
