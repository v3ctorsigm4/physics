
#include<stdio.h>
#include<stdlib.h>
#include<math.h>



int main()
{
	double mw=80.385,mz=91.1876,mh=125.09,iradian=M_PI/180.0,angw,e=1.6e-19,g,h,v,q,L;
mw*=e; // eV/c^2
mz*=e;

	angw=acos(mw/mz);
	angw=30.0*iradian;
	
	printf("%f \n",angw/iradian);
	g=sin(angw);	printf("%e g\n",g);//forditva!!!!!!! bug
	h=cos(angw);	printf("%e h\n",h);
printf("%e h+g\n",h+g);
printf("%e tan\n",atan(g/h)/iradian);
printf("%e cos\n",acos(g/h)/iradian);
	g=e/sin(angw);	printf("%e g\n",g);
	h=e/cos(angw);	printf("%e h\n",h);
//	printf("%e \n",atan(h/g)/iradian);


// MW = 1/2 vg
 v=2*mw/g; printf("%e v?\n",v);// 75.9??  v=246.0;//GeV
 L=mh*mh/(2.0*v*v); printf("%e L\n",L);// MH = sqrt(2Lv^2)
 q=atan(1.0/L)/iradian; printf("%e q\n",q);
 q=mh/sqrt(2.0*L*v*v);	printf("%e mh/mh\n",q);

 
	q=e/(g*h/sqrt(g*g+h*h));printf("%e e/e\n",q);
	
#if 0	
------------------------------	
a matrix	
	[2][0]
	[0][1]
eigenvalue
(A-IL)v=0

determinant =0
[a,b]
[c,d] = ad - bc =0	
	[2-L][0]
	[0][1-L]
 (2-L)*(1-L) -0 = 0
 (2-L -L(2-L) = 0
 (2 -3L + L^2) = 0
 
a=1
b=3
c=2
L=(-b+-sqrt(b^2 - 4ac)) /2a
L=(-3+-sqrt(3^2 - 42)) /2
L=(-3+-sqrt(9 - 8)) /2
L=(-3+-1) /2
eigenvalue
L1=-1
L2=-2

eigenvector
Av=Lv 
[2][0][v1] = L1[v1]
[0][1][v2]   L1[v2]
	
2*v1+0*v2 = L1*v1
0*v1+1*v2 = L1*v2
	
2*v1+0*v2 = L1*v1
0*v1+1*v2 = L1*v2

v1=1
2+0*v2 = L1
v2 = (L1-2/0)

v2=1
2*v1+0*1 = L1*v1
2 = L1 lol

------------------------------	
a matrix	
	[2][1]
	[1][2]
eigenvalue
(A-IL)v=0

determinant =0
[a,b]
[c,d] = ad - bc =0	
	[2-L][1]
	[1][2-L]
 (2-L)*(2-L) -1 = 0
 (4-2L)+ (2L-LL) -1 = 0
 3 -LL = 0
 
a=-1
b=0
c=3
L=(-b+-sqrt(b^2 - 4ac)) /2a
L=(+-sqrt( 12)) /-2
eigenvalue
L1=-3.4641 /2 = 1.732  = sqrt(3)
L2=3.4641 /2
{1/1.732=0.577367206 ?}

eigenvector
Av=Lv 
[2][1][v1] = L1[v1]
[1][2][v2]   L1[v2]
	
2*v1+1*v2 = L1*v1
1*v1+2*v2 = L1*v2
	
2*v1+v2 = L1*v1
 v1+2*v2 = L1*v2

v1=1
v2 = L1-2 = -1.732-2 = -3.732
check
2*1+1*-3.732 = L1*1
#endif


	return 0;
}











