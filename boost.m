close all 
clear all
L=50*10^(-6);
rL=10*10^(-3);
C=1000*10^(-6);
R=6.33;
rC=20*10^(-3);
Vin_max=15;
Vin_min=11.5;


alpha=R/(R+rC)


D=0:0.001:1;
res=[];
for i=1:1001
re=(R*(1-D(i)))/(alpha*rC*(1-D(i))+alpha*R*(1-D(i)).*(1-D(i))+rL);
res=[res re];
end
figure()
plot(D,res,'-');
% Vin = 15 => Duty cylce=21.3%
% Vin=11.5 => Duty cylce=40%

 

%gain max = 12.06 avec duty cycle D de 0.956
%sans prendre en compte les pertes 
D_max=1-(15/19)
D_min=1-(11.5/19)

% FT sans perte 
numerator=[-L*Vin_max R*(1-D_max)^2*Vin_max]
denominator=[(1-D_max)^2*C*L*R (1-D_max)^2*L R*(1-D_max)^4]
denominator2=[(1-D_min)^2*C*L*R (1-D_min)^2*L R*(1-D_min)^4]
numerator2=[-L*Vin_min R*(1-D_min)^2*Vin_min]

sys = tf(numerator,denominator)
%[A,B,C,D]=tf2ss(numerator,denominator);

sys2 = tf(numerator2,denominator2)
%[A2,B2,C2,D2]=tf2ss(numerator2,denominator2);

% plot des bodes des deux fonctions transfert 

 figure()
 bode(sys) 
 hold on 
 bode(sys2)


% Système d'états avec les pertes 
%Matrice et système perte Vmax
a00=rL+(1-D_max)*alpha*rC;
a01=(1-D_max)*alpha/L;
a10=((1-D_max)*alpha)/C;
a11=-alpha/(R*C);

A=[-a00/L -a01 ; 
   a10 a11] ;

V_bar1=alpha*rC*(1-D_max)+alpha*R*(1-D_max)*(1-D_max)+rL;

b00=Vin_max/(L*(1-D_max))
b01=-Vin_max/((1-D_max)*(1-D_max)*C*R);
B=[b00;
   b01];

C1=[alpha*rC*(1-D_max) alpha];


D_boost1=-V_bar1*(alpha*rC)/(R*(1-D_max));

%Matrice Perte Vmin

a00_min=rL+(1-D_min)*alpha*rC;
a01_min=(1-D_min)*alpha/L;
a10_min=((1-D_min)*alpha)/C;
a11_min=-alpha/(R*C);

A_min=[-a00_min/L -a01_min ; 
   a10_min a11_min] ;

V_bar_min=alpha*rC*(1-D_min)+alpha*R*(1-D_min)*(1-D_min)+rL;

b00_min=Vin_min/(L*(1-D_min))
b01_min=-Vin_min/((1-D_min)*(1-D_min)*C*R);
B_min=[b00_min;
   b01_min];

C1_min=[alpha*rC*(1-D_min) alpha];


D_boost_min=-V_bar_min*(alpha*rC)/(R*(1-D_min));



sys3=ss(A,B,C1,D_boost1)
hold on 
bode(sys3)


sys4=ss(A_min,B_min,C1_min,D_boost_min)
hold on 
bode(sys4)



 legend("sans perte V_in=15V","sans perte V_in=11.5V","avec perte V_in=15V","avec perte V_in=11.5V")
 
 
 C1=[alpha*rC alpha];
 C2=[0 alpha];
 A1=[-(rL+alpha*rC)/L -alpha/L;
     alpha/C -alpha/(R*C)];
 A2=[-rL/L 0;
     0 -alpha/(R*C)];
 B1=[1/L;
     0];
 
 
 % Paramètre correcteur avance de phase 
 
[Num,Denum]=ss2tf(A,B,C1,D_boost1)

Denum2=conv([1 0],Denum)
sys_integre=tf(Num,Denum2)
figure()
bode(sys_integre)
% frequence entre 10 200 rad/sec et 75 000rad/sec => 20 000 rad/sec 
[Gm,Pm,Wcg,Wcp] = margin(sys_integre)

% marge de phase  89.9

K=7; 
 
 phiA=60-(-90);
 a=(1+sin(phiA/180*pi))/(1-sin(phiA/180*pi))
% calcul wm fréquence quand le gain est égal à -10log(a)
wm=2e4
%T1 et T2 égaux à T
T1=1/(wm*sqrt(a))
T2=1/(wm*sqrt(a))