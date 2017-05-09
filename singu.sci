clear
function c=singu0(x,a)
    c=x>=a;
endfunction

function c=singu1(x,a)
    c=(x-a).*(x>=a);
endfunction

function c=singu2(x,a)
    c=((x-a).^2).*(x>=a);
endfunction

function c=singu3(x,a)
    c=((x-a).^3).*(x>=a);
endfunction

function c=singu4(x,a)
    c=((x-a).^4).*(x>=a);
endfunction

function dibujar()
    subplot(2,2,1)
    plot2d(x,Vx);
    a=gca(); // Handle on axes entity
    a.x_location = "origin";
    a.y_location = "origin";
    e = gce();
    e.children.thickness = 3;
    xtitle("Fuerza cortante");
    xlabel("x")
    ylabel("V(x)")
    subplot(2,2,3)
    plot(x,Mx);
    a=gca(); // Handle on axes entity
    a.x_location = "origin";
    a.y_location = "origin";
    e = gce();
    e.children.thickness = 3;
    p=e.children(1);
    [m,k]=max(Mx)
    t=datatipCreate(p,k);
    xtitle("Momento flector");
    xlabel("x");
    ylabel("M(x)")
    subplot(2,2,2);
    plot(x,theta);
    a=gca(); // Handle on axes entity
    a.x_location = "origin";
    a.y_location = "origin";
    e = gce();
    e.children.thickness = 3;
    xtitle("Angulo de pendiente");
    xlabel("x");
    ylabel("Theta(x)")
    subplot(2,2,4);
    plot(x,y);
    a=gca(); // Handle on axes entity
    a.x_location = "origin";
    a.y_location = "origin";
    e = gce();
    e.children.thickness = 3;
    xtitle("Deflexion");
    xlabel("x");
    ylabel("y(x)")
endfunction

//PARAMETROS DE ENTRADA
L=10;               //metros
F=(91/3)*10^3;      //Newtons
W=24375;            //Newtons por metro
M=2.75;             //metros
S=3.8;              //metros
I=462*10^-6;        //metros^4
E=200*10^9;         //Pascales
K1=%inf;          //Newtons/metro
K2=%inf;          //Newtons/metro     %inf
d=420*10^(-3);      //metros
sadm=250*10^6;      //Pascales
//PARAMETROS DE ENTRADA

x=0:0.001:L;
   
   if K1==0 then
       C1=1;
   else
       C1=K1;
   end
   
   if K2==0 then
       C2=1;
   else
       C2=K2;
   end
   
A=[1,0,1,1;
   0,1,L,S;
   (S^3)/6,-(S^2)/2,0,((E*I)/C1);
   (L^3)/6,-(L^2)/2,((E*I)/C2),((L-S)^3)/6];
   
   

b=[-F-W*(L-S);
    F*(-1-M)-W*(L-S)*(((L-S)/2)+S);
    F*(((S-M)^2/2)-(S-M)^3/6);
    F*(((L-M)^2/2)-(L-M)^3/6)-W*(L-S)^4/24];
    
    if K1==0 then
       A(3,:)=0
       A(:,4)=0
       b(3)=0
   end
   
   if K2==0 then
       A(4,:)=0
       A(:,3)=0
       b(4)=0
   end

[x0,nsA]=linsolve(A,b);

Fa=x0(1,1);
Fc=x0(4,1);
Fd=x0(3,1);
Mr=x0(2,1);

Vx=Fa*singu0(x,0)-F*singu0(x,M)+Fc*singu0(x,S)-W*singu1(x,S)+W*singu1(x,L)+Fd*singu0(x,L);
Mx=Fa*singu1(x,0)-F*singu1(x,M)+Fc*singu1(x,S)-(W/2)*singu2(x,S)+(W/2)*singu2(x,L)+Fd*singu1(x,L)-Mr*singu0(x,0)+F*singu0(x,M);

sigma=(max(abs(Mx))*d/2)/I;
Fs=sadm/sigma;

function theta=angulo(x)
    theta=((Fa/2)*singu2(x,0)-(F/2)*singu2(x,M)+(Fc/2)*singu2(x,S)-(W/6)*singu3(x,S)+(W/6)*singu3(x,L)+(Fd/2)*singu2(x,L)-Mr*singu1(x,0)+F*singu1(x,M))/(E*I);
endfunction
theta=angulo(x);

function y=deflex(x)
    y=((Fa/6)*singu3(x,0)-(F/6)*singu3(x,M)+(Fc/6)*singu3(x,S)-(W/24)*singu4(x,S)+(W/24)*singu4(x,L)+(Fd/6)*singu3(x,L)-(Mr/2)*singu2(x,0)+(F/2)*singu2(x,M))/(E*I);
endfunction
y=deflex(x);

resp=zeros(2,3);

resp(1,1)=angulo(M); resp(1,2)=angulo(S); resp(1,3)=angulo(L);
resp(2,1)=deflex(M); resp(2,2)=deflex(S); resp(2,3)=deflex(L);

mprintf('\n Para las variables de entrada el analisis resultante es el siguiente:')
mprintf('\n')
mprintf('\n Las reacciones en los apoyos son: \n Fa=%1.3f kN \n Fc=%1.3f kN \n Fd=%1.3f kN \n Mr=%1.3f kN*m',Fa/1000,Fc/1000,Fd/1000,Mr/1000)

if sigma>sadm then
    mprintf('\n')
    mprintf('\n -_-_-_-_-_-¡ALERTA!-_-_-_-_-_-')
    mprintf('\n Se ha sobrepasado el esfuerzo normal admisible en la viga')
    mprintf('\n Esfuerzo admisible [Õ]=%1.3f MPa',sadm/1000000)
    mprintf('\n Esfuerzo maximo en la viga Õ=%1.3f MPa',sigma/1000000)
end

mprintf('\n')
mprintf('\n Angulo de pendiente en el punto B=%1.2f metros es Theta=%f radianes',M,resp(1,1))
mprintf('\n Angulo de pendiente en el punto C=%1.2f metros es Theta=%f radianes',S,resp(1,2))
mprintf('\n Angulo de pendiente en el punto D=%1.2f metros es Theta=%f radianes',L,resp(1,3))

mprintf('\n')
mprintf('\n Deflexión en el punto B=%1.2f metros es y=%1.3f mm',M,(resp(2,1)*1000))
mprintf('\n Deflexión en el punto C=%1.2f metros es y=%1.3f mm',S,(resp(2,2)*1000))
mprintf('\n Deflexión en el punto D=%1.2f metros es y=%1.3f mm',L,(resp(2,3)*1000))

mprintf('\n')
mprintf('\n Con un momento máximo de M=%1.3f kN*m y con la caracteristicas \n de la seccion de la viga se presenta un esfuerzo \n normal maximo de Õ=%1.3f MPa',max(abs(Mx))/1000,sigma/1000000 )
mprintf('\n De esta manera la viga presenta un factor de seguridad Fs=%1.3f',Fs)
