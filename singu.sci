clear
close();
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
    t=a.title;
    t.foreground=9;
    t.font_size=4;
    t.font_style=6;
    t=a.y_label;
    t.font_size=3;
    //momento flector
    subplot(2,2,3)
    plot(x,Mx);
    a=gca(); // Handle on axes entity
    a.x_location = "origin";
    a.y_location = "origin";
    e = gce();
    e.children.thickness = 3;
    p=e.children(1);
    [m,k]=max(abs(Mx))
    t=datatipCreate(p,k);
    xtitle("Momento flector");
    xlabel("x");
    ylabel("M(x)")
    t=a.title;
    t.foreground=9;
    t.font_size=4;
    t.font_style=6;
    t=a.y_label;
    t.font_size=3;
    //angulo deflexion
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
    t=a.title;
    t.foreground=9;
    t.font_size=4;
    t.font_style=6;
    t=a.y_label;
    t.font_size=3;
    //deflexion
    plot(x,y);
    a=gca(); // Handle on axes entity
    a.x_location = "origin";
    a.y_location = "origin";
    e = gce();
    e.children.thickness = 3;
    xtitle("Deflexion");
    xlabel("x");
    ylabel("y(x)")
    t=a.title;
    t.foreground=9;
    t.font_size=4;
    t.font_style=6;
    t=a.y_label;
    t.font_size=3;
endfunction

//PARAMETROS DE ENTRADA
L=10;               //metros
F=(91/3)*10^3;      //Newtons
W=24375;            //Newtons por metro
M=2.75;             //metros
S=3.8;              //metros
    //Parametros de sección especificos
    tf=19.3*10^(-3);
    bf=261*10^(-3);
    d=420*10^(-3);
    tw=11.6*10^(-3);
    I=2*((1/12*bf*tf^3)+bf*tf*(d/2-tf/2)^2)+(1/12*tw*(d-2*tf)^3)
    //Parametros de sección especificos
//I=462*10^-6;        //metros^4
E=200*10^9;         //Pascales
K1=0;          //Newtons/metro     %inf para valor infinito
K2=0;          //Newtons/metro     %inf para valor infinito
sadm=250*10^6;      //Pascales
thetadm=3;          //grados
yadm=10*10^(-3);    //milimetros
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

function theta=angulo(x)
    theta=(((Fa/2)*singu2(x,0)-(F/2)*singu2(x,M)+(Fc/2)*singu2(x,S)-(W/6)*singu3(x,S)+(W/6)*singu3(x,L)+(Fd/2)*singu2(x,L)-Mr*singu1(x,0)+F*singu1(x,M))/(E*I))*180/%pi;
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

sigmax=(max(abs(Mx))*d/2)/I;
Fs=sadm/sigmax;
if sigmax>sadm then
    mprintf('\n')
    mprintf('\n -_-_-_-_-_-¡ALERTA!-_-_-_-_-_-')
    mprintf('\n Se ha sobrepasado el esfuerzo normal admisible en la viga')
    mprintf('\n Esfuerzo admisible [Õ]=%1.3f MPa',sadm/1000000)
    mprintf('\n Esfuerzo maximo en la viga Õ=%1.3f MPa',sigmax/1000000)
end

[thetamax,k]=max(abs(theta));
if thetamax>thetadm then
    mprintf('\n')
    mprintf('\n -_-_-_-_-_-¡ALERTA!-_-_-_-_-_-')
    mprintf('\n Se ha sobrepasado el angulo de deflexión admisible en la viga')
    mprintf('\n Angulo admisible [Theta]=%f Grados',thetadm)
    mprintf('\n Angulo maximo de la viga Theta=%f grados que se da en x=%1.2f metros',thetamax,(k-1)*0.001)
end

[ymax,k]=max(abs(y));
if ymax>yadm then
    mprintf('\n')
    mprintf('\n -_-_-_-_-_-¡ALERTA!-_-_-_-_-_-')
    mprintf('\n Se ha sobrepasado la deflexión admisible en la viga')
    mprintf('\n Deflexión admisible [y]=%1.3f mm',yadm*1000)
    mprintf('\n La deflexion maxima de la viga es y=%f mm y se da en x=%1.2f metros',ymax*1000,(k-1)*0.001)
end

mprintf('\n')
mprintf('\n Angulo de pendiente en el punto B=%1.2f metros es Theta=%f grados',M,(resp(1,1)))
mprintf('\n Angulo de pendiente en el punto C=%1.2f metros es Theta=%f grados',S,resp(1,2))
mprintf('\n Angulo de pendiente en el punto D=%1.2f metros es Theta=%f grados',L,resp(1,3))

mprintf('\n')
mprintf('\n Deflexión en el punto B=%1.2f metros es y=%1.3f mm',M,(resp(2,1)*1000))
mprintf('\n Deflexión en el punto C=%1.2f metros es y=%1.3f mm',S,(resp(2,2)*1000))
mprintf('\n Deflexión en el punto D=%1.2f metros es y=%1.3f mm',L,(resp(2,3)*1000))

mprintf('\n')
mprintf('\n Con un momento máximo de M=%1.3f kN*m y con la caracteristicas \n de la seccion de la viga se presenta un esfuerzo \n normal maximo de Õ=%1.3f MPa',max(abs(Mx))/1000,sigmax/1000000 )
mprintf('\n De esta manera la viga presenta un factor de seguridad Fs=%1.3f',Fs)
 dibujar;
