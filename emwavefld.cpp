 


#include "eng6.cpp"


//#define _EONLY
#define _EandB
//#define _SKYRME2

// NA
//#define _BONLY
//#define _AXION
//#define _SKYRME


//#define _TESTCHARGE
//#define _EDUMP



#define DXX 0.01
//0.01

/*
BUGS
a t=1234 nem jo mert trafo utal elmaszik
ha cheatx and _INVMAT akkor jo minden elojel!!!,    / mar nincs/

Tdir transformalt is kell
Tdir_tr  , de hossz kontrakcio nelkul
time mindig kell
A vektorpotencial inverzben transzformalodik         (qmatInverse =lorentz_transformation2i)
skalarmezo minusz
dx es x ugyanabban az IR-ben van!, szoval egyszerre 
-gradFI +dA
curl +
/dx*2.0 mindig


transformEMfield elotrafo nelkul is rossz a szog
*/


// WM wave field
// A vector pot wave csak akkor jo ha X componens van  LOL  ,    de a tranformalt dt derivalt valtozik!!
// S scalar wave     fncEgrad() dx2 vel JO!!!!
// static = emcurl2c.cpp      EMwave in X emwavefld.cpp


float1 k=M_PI/180.0;//80
float1 c=1.0;//3e8; ua


float1 v_def=0.85;
float1 v=v_def*c;//0.3-0.9
//float1 v=0.0*c;
float1 gamma2, beta,time2=0.0;//       gamma used in gcc!!!!!!

//FREE       barmilyen iranyjo!!!!      xyz jo kulon, keveredve kicsit elfordul
//#ifdef _EONLY
#if 0
//3: 77,72 fok  90 OK
//vec3 move_dir=normalize(vec3(0,10,20));// ez iranyu komponens az X   0,10,20
//vec3 Tdir=normalize(vec3(1.0,1.0,0.0));	//nem Xbe megy       hullam,wave
//vec3 move_dir=normalize(vec3(2,15,-11));vec3 Tdir=normalize(vec3(1.0,0.3,0.6));	
vec3 move_dir=normalize(vec3(10,0,0));vec3 Tdir=normalize(vec3(1.0,1.0,0.0));	
#else
//3: 86,85 fok  90 OK
vec3 move_dir=normalize(vec3(10,0,0));//45 fok
vec3 Tdir=normalize(vec3(1.0,0.0,1.0));	
#endif
vec3 Tdir_tr,move_dir_tr;


qmat tensorF;//  [4][4];// EM tensor    , The electromagnetic field strength tensor
qmat tensorV;// jele:V fejjel!!      boost matrix
qmat tensorF2;// F' = transformed F
qmat tensorF3;//temp



int espacetime=0;

//--------------------------------------------------------------------------------------

vec3 lorentz_transformation2(float1 t1,vec3 P1,float1 &t2,int opt=0)//t1 == time component of fourvector
{
	vec3 P2=P1;//	float1 gamma2=1.0/sqrt(1.0-v*v/(c*c));
	float1 sg=1.0,x2,x1=dot(P1,move_dir);//      move_dir_tr UA!
	if(opt&4) t1=-t1;// -   S mindig minus lorentz elott     OK   -+++
	
//USED == 0,1,2,16	
	if(opt&16) sg=-1.0;//inv transf   0x10
	
	if(opt&1)//no length contraction  ,           (csak a Tdit_tr -nel !)
	{
		x2=(x1 - sg*v*t1);// no length contr
		t2=(t1 - sg*v*x1/(c*c));
	}
	else//0-2
	{
		x2=(x1 - sg*v*t1)*gamma2;//Lorentz transformation                   mashol igy van
		t2=(t1 - sg*v*x1/(c*c))*gamma2;  //       a transzformalt skalar field
	}

	//extension?//NA	
	if(opt&2) x2-=sg*v*t2;// we need t2==0 position              ezzel meg rosszabb

	P2 -= move_dir*x1;
	P2 += move_dir*x2;

	if(opt&4) t2=-t2;//! -+++

	return P2;
}


vec3 lorentz_transformation4(float1 t1,vec3 P1,float1 &t2,int inv=0);
//UA ugyanaz a letto OK!        (ha -time)            2 vs 4
#define lorentz_transformationA lorentz_transformation4
#define lorentz_transformationB lorentz_transformation4


int LT1=16;//16 0 OLD
int LT2=0;//inv(16) + sign(4)



//--------------------------------------------------------------------------------------


float1 fncFI_base(float1 x,float1 y,float1 z,float1 t,int tra=0)//scalar potential  =  FI
{
	vec3 Tdir3=Tdir; if(tra) Tdir3=Tdir_tr;

	vec3 Pos=vec3(x,y,z);                           
	float1 dist=dot(Pos,Tdir3) -t;//length(P)    //sose _tr ?        +time!

	float1 phase=dist*k;
#ifdef _EONLY
//	vec3 Pos2=vec3(0,y,z);  	return -1e3/length(Pos2);//scalar field            current?
	return -1e3/length(Pos);//scalar field            static
//	return 1e6/dot(Pos,Pos);//scalar field            static
//	return (1e3+100.0*t)/length(Pos);//scalar field   dynamics
#endif
#ifdef _EandB
//	return cos(phase)*0.5;// origin
	//return cos(phase)*0.965925826;//15
	//return 0.0;
	
	float1 rotang=65.0*iradian;	
	return cos(phase);//*cos(rotang);// OK!     sin/cos fuggoen amilyen a Tdir    -??
#endif
#ifdef _BONLY
	return 0.0;//none
#endif	
#ifdef _AXION
	return 1.0;// vortex, orveny,   magnetic monopole? nem LOL
#endif	
#ifdef _SKYRME
	{
	vec3 axis=normalize(vec3(x,0.0,z))*150.0;// torus center ring
	vec3 side=normalize(cross(axis,vec3(0,1,0))); // UP
side+=normalize(axis); side=normalize(side);	
//	vec3 dPos=vec3(x,t,z)-axis; //  t==y
	vec3 dPos=vec3(x,y,z)-axis; // igy toltes
//	float1 dist=length4(dPos,y) ;
	float1 dist=length(dPos) ;
	vec3 A=normalize(cross(dPos,side))*20000.0/(dist*dist);
	return A.y;// y=t
	}
#endif	
#ifdef _SKYRME2
	return 0.0;// none
#endif	


//	return 1e0/(l);// elvileg 1/x^3
}



vec3 fncA_base(float1 x,float1 y,float1 z,float1 t,int tra=0)// A4 = (FI/c, A.x, A.y, A.z)
{
	vec3 Tdir3=Tdir; if(tra) Tdir3=Tdir_tr;
	vec3 A;

#ifdef _AXION
	{
	vec3 Pos=vec3(x,y,z);// AXION  magnetic
	float1 dist=length(Pos) ;//+t;
	A=normalize(vec3(-z,0.0,x))*150000.0/(dist*dist);// vortex, orveny,   magnetic monopole? dipole LOL
	return A;
	}
#endif

#ifdef _SKYRME
	{//4d     van E   de nincs B  xD
	vec3 axis=normalize(vec3(x,0.0,z))*150.0;// torus center ring
	vec3 side=normalize(cross(axis,vec3(0,1,0))); // UP
side+=normalize(axis); side=normalize(side);
//	vec3 dPos=vec3(x,t,z)-axis; // y==t
	vec3 dPos=vec3(x,y,z)-axis;  // igy +sima length jo, de logikatlan
//	float1 dist=length4(dPos,y) ;//
	float1 dist=length(dPos) ;//
	A=normalize(cross(dPos,side))*20000.0/(dist*dist);
	A.y=0.0;// into time dimension        t=y lesz!
	return A;
	}
#endif
#ifdef _SKYRME2
	{//sima 3d , nincs E ter!
	vec3 axis=normalize(vec3(x,0.0,z))*150.0;// torus center ring
	vec3 side=normalize(cross(axis,vec3(0,1,0))); // UP
	vec3 dPos=vec3(x,y,z)-axis; 
	float1 dist=length(dPos) ;//+t;
	A=normalize(cross(dPos,side))*140000.0/(dist*dist);//20000  /t^2
	return A;
	}
#endif

	
	
#ifdef _EONLY
{
	A=vec3(0,0,0);//  Csak scalar wave
vec3 Pos=vec3(x,y,z);
float1 dist=length(Pos);
vec3 norm=normalize(cross(Pos,vec3(0,1,0)));//y
A=norm*90000.0/(dist*dist);

	return A;
}
#endif

//#ifdef _BONLY
	vec3 Pos=vec3(x,y,z);//
	float1 dist=dot(Pos,Tdir3) -t;//delta == time derivative       //sose _tr ?    +time!
	float1 phase=dist*k;
#if 1
//Tdir3-nal csak ez lehet !!! mert elfordul!
	vec3 Tdir2; Tdir2.x=Tdir3.y;Tdir2.y=Tdir3.z;Tdir2.z=Tdir3.x;//  ferde
	vec3 Xdir=cross(Tdir3,Tdir2);
	vec3 Ydir=normalize(cross(Xdir,Tdir3));
	Xdir =normalize(cross(Tdir3,Ydir));
#else//most nem lehet X mert elol van Tdir!!!	
//	vec3 Xdir=normalize(vec3(1.0, 0.0,0.0));// x-be jo  	
	vec3 Xdir=normalize(vec3(0.0, 0.0,1.0));// z-be nem!!!!!!!!11  mivel nincs xkomponens LOOL
#endif

	//Xdir=normalize(vec3(0.0, 0.0,1.0));// z-be nem!!!!!!!!11  mivel nincs xkomponens LOOL
	A=Xdir*sin(phase);//mindig kell, mert korbe kell mennie,   (a derivalas miatt)
	
#ifdef _EandB
//barmilyen ellenvektor jo, ha a mozgasiranyu komponense ellentetes (most azonos, mivel elojelhiba van)
//	float1 rotang=acos(0.5/0.866);//  ferde
	float1 rotang=65.0*iradian;
//	A+=Tdir3*cos(phase)*cos(rotang);       //    A()=Tdir3*cos(phase)*cos(rotang)*0.866 amp=0.5 , S()=-cos(phase)*0.5 , amp=0.5
//	A+=Ydir*cos(phase)*sin(rotang);// akkor forog, ha ez van
	A+=Tdir3*cos(phase);       //    A()=Tdir3*cos(phase)*cos(rotang)*0.866 amp=0.5 , S()=-cos(phase)*0.5 , amp=0.5
//	A+=Ydir*cos(phase);//??


//A*=(1.0+t*0.1);
	
/*	A+=Tdir3*cos(phase)*cos(rotang);// origin
	A+=Ydir*cos(phase)*sin(rotang);
	A*=0.866;*/
#endif	
//#endif

	return A;
}





vec3 fncA_tr(float1 x,float1 y,float1 z,float1 t)
{
	float1 t2=0.0;//                                length contraction
vec3 P(x,y,z);vec3 P2=lorentz_transformationA(t,P,t2,LT1);x=P2.x;y=P2.y;z=P2.z; t=t2;  //field position!  2?3

	float1 S=fncFI_base(x,y,z,t,1);//FI idofuggo!        mimdig kell time!!!!!!! 
	vec3 A=fncA_base(x,y,z,t,1);
	vec3 A2=lorentz_transformationB(S/c,A,t2,LT2);  

	return A2;
}
float1 fncFI_tr(float1 x,float1 y,float1 z,float1 t)
{
	float1 t2=0.0;//                                length contraction
vec3 P(x,y,z);vec3 P2=lorentz_transformationA(t,P,t2,LT1);x=P2.x;y=P2.y;z=P2.z; t=t2;  //field position!  2?3

	float1 S=fncFI_base(x,y,z,t,1);// ha mar transformalt, akkor Tdir_tr kell!!!   
	vec3 A=fncA_base(x,y,z,t,1);
	vec3 A2=lorentz_transformationB(S/c,A,t2,LT2);  // inverse matrix
	
	return t2*c;
}
//press 4
vec4 fncA4_base(float1 t,float1 x,float1 y,float1 z)//  !!!!!!!! txyz
{
// 1kell mivel fncdAtransfdA mar transformal!   , DE a transformEMfield sem Tdir_tr hasznal!!!

	float1 S=fncFI_base(x,y,z,t,0);// 0 OK     
	vec3 A2=fncA_base(x,y,z,t,0);

	vec4 A4;
	A4.x=A2.x;
	A4.y=A2.y;
	A4.z=A2.z;
	A4.w=S; //    FI idofuggo!
	
	return A4;
}
//press 5 OK
vec4 fncA4_tr(float1 t,float1 x,float1 y,float1 z)//  !!!!!!!! txyz
{
	float1 t2=0.0;//                                  length contraction
if(v>0.0)
{
	vec3 P(x,y,z);vec3 P2=lorentz_transformationA(t,P,t2,LT1);x=P2.x;y=P2.y;z=P2.z; t=t2;  //field position!  2?3
}
else t2=t;
	float1 S=fncFI_base(x,y,z,t,1);// 1=masik RF ben vagyunk! Tdir_tr 
	vec3 A=fncA_base(x,y,z,t,1);
	vec3 A2=lorentz_transformationB(S/c,A,t2,LT2); 
if(v==0.0) A2=A;
		
	vec4 A4;
	A4.x=A2.x;
	A4.y=A2.y;
	A4.z=A2.z;
	A4.w=t2; 
	
	return A4;
}

//--------------------------------------------------------------------------------------
//ok
vec3 fncEgrad_base(float1 x,float1 y,float1 z,float1 t)//E=gradient of scalar potential    E = - grad FI
{
	float1 dx=DXX, dy=DXX, dz=DXX ,dt=DXX;//  dx es x ugyanabban az IR-ben van!
	
	float1 dSx=(fncFI_base(x+dx,y   ,z	 ,t   ) - fncFI_base(x-dx,y   ,z   ,t   ))/(dx*2.0);	// OK
	float1 dSy=(fncFI_base(x   ,y+dy,z	 ,t   ) - fncFI_base(x   ,y-dy,z   ,t   ))/(dy*2.0);
	float1 dSz=(fncFI_base(x   ,y   ,z+dz,t   ) - fncFI_base(x   ,y   ,z-dz,t   ))/(dz*2.0);
	vec3    dA=( fncA_base(x   ,y   ,z   ,t+dt) -  fncA_base(x   ,y   ,z   ,t-dt))/(dt*2.0); // +-!  cheatx

	vec3 gradFI=vec3(dSx, dSy, dSz); //OK
	
	//same as four gradient:  E=idoderival(dA) - scalar derivaltja(gradFI)
	return -gradFI -dA;//  gradient of FI      CHEAT +dA       --OK
}

vec3 fncEgrad_tr(float1 x,float1 y,float1 z,float1 t)//E=gradient of scalar potential    E = - grad FI
{
	float1 dx=DXX, dy=DXX, dz=DXX, dt=DXX;
	
	float1 dSx=(fncFI_tr(x+dx,y   ,z   ,t   ) - fncFI_tr(x-dx,y   ,z   ,t   ))/(dx*2.0);	// OK
	float1 dSy=(fncFI_tr(x   ,y+dy,z   ,t   ) - fncFI_tr(x   ,y-dy,z   ,t   ))/(dy*2.0);
	float1 dSz=(fncFI_tr(x   ,y   ,z+dz,t   ) - fncFI_tr(x   ,y   ,z-dz,t   ))/(dz*2.0);
	vec3    dA=( fncA_tr(x   ,y   ,z   ,t+dt) -  fncA_tr(x   ,y   ,z   ,t-dt))/(dt*2.0); //  OK  +-!  cheatx

	vec3 gradFI=vec3(dSx, dSy, dSz); //OK
	
	return -gradFI -dA;//  gradient of FI      CHEAT +dA !!!  --OK
}
//pedig E.x=dAw.x - dAx.w;
//szoval a time derivalt x je (dA.x) + 
//a scalarfield(w) x derivaltja (gradFI.x=Sdx-Sdx2 szoval elvileg jo)



vec3 curlA_tr(float1 x,float1 y,float1 z,float1 t)	
{
	float1 dx=DXX, dy=DXX, dz=DXX;

	vec3 dAx=(fncA_tr(x+dx,y,z,t) - fncA_tr(x-dx,y,z,t))/(dx*2.0);  // ok
	vec3 dAy=(fncA_tr(x,y+dy,z,t) - fncA_tr(x,y-dy,z,t))/(dy*2.0);
	vec3 dAz=(fncA_tr(x,y,z+dz,t) - fncA_tr(x,y,z-dz,t))/(dz*2.0);

//B= curl A
	vec3 curl= vec3((dAz.y)-(dAy.z),//def X f         -?    nem minusznak kell lennie!       CHEAT  +OK
					(dAx.z)-(dAz.x),	//like cross product
					(dAy.x)-(dAx.y));
					
	return curl;
}
vec3 curlA_base(float1 x,float1 y,float1 z,float1 t)	
{
	float1 dx=DXX, dy=DXX, dz=DXX;

	vec3 dAx=(fncA_base(x+dx,y,z,t) - fncA_base(x-dx,y,z,t))/(dx*2.0);// ok
	vec3 dAy=(fncA_base(x,y+dy,z,t) - fncA_base(x,y-dy,z,t))/(dy*2.0);
	vec3 dAz=(fncA_base(x,y,z+dz,t) - fncA_base(x,y,z-dz,t))/(dz*2.0);

//B= curl A
	vec3 curl= vec3((dAz.y)-(dAy.z),//def X f      CHEAT  
					(dAx.z)-(dAz.x),	//like cross product
					(dAy.x)-(dAx.y));
					
	return curl;
}

//--------------------------------------------------------------------------------------


void tensorVVFProduct()//OKE!       F -> F2
{
//mat[a][b]=mat1[a][i]*mat2[b][j]*mat3[i][j]

	for(int a=0;a<4;a++)//256 step LOL
	for(int b=0;b<4;b++)
	{
		tensorF2.m[a][b]=0.0;
		
		for(int i=0;i<4;i++)
		for(int j=0;j<4;j++)
		{
			tensorF2.m[a][b] += (tensorV.m[a][i]*tensorV.m[b][j]*tensorF.m[i][j]);
//			tensorF2.m[a][b] += (tensorV.m[b][i]*tensorV.m[a][j]*tensorF.m[j][i]);//NA  UA! ( mivel V szimmetrikus)
		}
	}
}

//press 4
void buildTensorF_base(float1 x,float1 y,float1 z,float1 t)
{
	vec4 dAx,dAy,dAz,dAw;
	float1 dt=DXX, dx=DXX, dy=DXX, dz=DXX;// -dt ? NO!

//ha -dt akkor fncEgrad_base miert nem ugyanugy?

//four gradient  -+++
	dAw=(fncA4_base(t-dt,x   ,y   ,z   ) - fncA4_base(t+dt,x   ,y   ,z   ))/(dt*2.0);//-+++   itt - , fncA +  ok   +-! cheatx
	dAx=(fncA4_base(t   ,x+dx,y   ,z   ) - fncA4_base(t   ,x-dx,y   ,z   ))/(dx*2.0);
	dAy=(fncA4_base(t   ,x   ,y+dy,z   ) - fncA4_base(t   ,x   ,y-dy,z   ))/(dy*2.0);
	dAz=(fncA4_base(t   ,x   ,y   ,z+dz) - fncA4_base(t   ,x   ,y   ,z-dz))/(dz*2.0);
	
//E.z =0 ,     z in w =0        w in z ==0 because Tdir.z=0,   de dAz.w!=0 !  vec4-operator BUG!
//float1 q1=fncA4(w   ,x   ,y   ,z+dx).w;float1 q2=fncA4(w   ,x   ,y   ,z-dx).w;


	tensorF.m[0][0]= (dAw.w - dAw.w);// w == time
	tensorF.m[0][1]= (dAw.x - dAx.w);//idoderival - scalar derivaltja
	tensorF.m[0][2]= (dAw.y - dAy.w);
	tensorF.m[0][3]= (dAw.z - dAz.w);

	tensorF.m[1][0]= (dAx.w - dAw.x);
	tensorF.m[1][1]=-(dAx.x - dAx.x);// - because of Hodge star (wedge vs cross product)
	tensorF.m[1][2]=-(dAx.y - dAy.x);
	tensorF.m[1][3]=-(dAx.z - dAz.x);

	tensorF.m[2][0]= (dAy.w - dAw.y);
	tensorF.m[2][1]=-(dAy.x - dAx.y);//- Hodge*
	tensorF.m[2][2]=-(dAy.y - dAy.y);
	tensorF.m[2][3]=-(dAy.z - dAz.y);

	tensorF.m[3][0]= (dAz.w - dAw.z);
	tensorF.m[3][1]=-(dAz.x - dAx.z);//- Hodge*
	tensorF.m[3][2]=-(dAz.y - dAy.z);
	tensorF.m[3][3]=-(dAz.z - dAz.z);
#if 0
	printf("%.3f %.3f %.3f %.3f \n",dAw.w,dAw.x,dAw.y,dAw.z);
	printf("%.3f %.3f %.3f %.3f \n",dAx.w,dAx.x,dAx.y,dAx.z);
	printf("%.3f %.3f %.3f %.3f \n",dAy.w,dAy.x,dAy.y,dAy.z);
	printf("%.3f %.3f %.3f %.3f \n",dAz.w,dAz.x,dAz.y,dAz.z);
	printf("\n");
	printf("%.3f %.3f %.3f %.3f \n",tensorF.m[0][0],tensorF.m[0][1],tensorF.m[0][2],tensorF.m[0][3]);
	printf("%.3f %.3f %.3f %.3f \n",tensorF.m[1][0],tensorF.m[1][1],tensorF.m[1][2],tensorF.m[1][3]);
	printf("%.3f %.3f %.3f %.3f \n",tensorF.m[2][0],tensorF.m[2][1],tensorF.m[2][2],tensorF.m[2][3]);
	printf("%.3f %.3f %.3f %.3f \n",tensorF.m[3][0],tensorF.m[3][1],tensorF.m[3][2],tensorF.m[3][3]);
	printf("\n");
#endif
}
//press 5 OK
void buildTensorF_tr(float1 x,float1 y,float1 z,float1 t)//_Tr nel mindig atvitel van time-ba
{
	vec4 dAx,dAy,dAz,dAw;
	float1 dt=DXX, dx=DXX, dy=DXX, dz=DXX;// -dt ? NO!
//transformed already

//four gradient     -+++
	dAw=(fncA4_tr(t-dt,x   ,y   ,z   ) - fncA4_tr(t+dt,x   ,y   ,z   ))/(dt*2.0);//-+++   itt - , fncA +  ok   +-! cheatx
	dAx=(fncA4_tr(t   ,x+dx,y   ,z   ) - fncA4_tr(t   ,x-dx,y   ,z   ))/(dx*2.0);
	dAy=(fncA4_tr(t   ,x   ,y+dy,z   ) - fncA4_tr(t   ,x   ,y-dy,z   ))/(dy*2.0);
	dAz=(fncA4_tr(t   ,x   ,y   ,z+dz) - fncA4_tr(t   ,x   ,y   ,z-dz))/(dz*2.0); // 
	
//E.z =0 ,     z in w =0        w in z ==0 because Tdir.z=0,   de dAz.w!=0 !  vec4-operator BUG!
//float1 q1=fncA4(w   ,x   ,y   ,z+dx).w;float1 q2=fncA4(w   ,x   ,y   ,z-dx).w;


	tensorF.m[0][0]= (dAw.w - dAw.w);// w == time
	tensorF.m[0][1]= (dAw.x - dAx.w);//Ex   idoderival - scalar derivaltja
	tensorF.m[0][2]= (dAw.y - dAy.w);//Ey
	tensorF.m[0][3]= (dAw.z - dAz.w);//Ez

	tensorF.m[1][0]= (dAx.w - dAw.x);
	tensorF.m[1][1]=-(dAx.x - dAx.x);// - because of Hodge star (wedge vs cross product)
	tensorF.m[1][2]=-(dAx.y - dAy.x);//Bz
	tensorF.m[1][3]=-(dAx.z - dAz.x);

	tensorF.m[2][0]= (dAy.w - dAw.y);
	tensorF.m[2][1]=-(dAy.x - dAx.y);
	tensorF.m[2][2]=-(dAy.y - dAy.y);
	tensorF.m[2][3]=-(dAy.z - dAz.y);//Bx

	tensorF.m[3][0]= (dAz.w - dAw.z);
	tensorF.m[3][1]=-(dAz.x - dAx.z);//By  elvileg ([0][0] [1][2] [2][3] [3][1] egy 4d CURL , es 4 van belole)
	tensorF.m[3][2]=-(dAz.y - dAy.z);
	tensorF.m[3][3]=-(dAz.z - dAz.z);

}
#if 0
//NA
	for(int i=0;i<4;i++)
	for(int j=0;j<4;j++) tensorF.m[i][j]=-tensorF.m[i][j];//jo ha mashol(fncFI) nincs -
#endif	



void fncdA(float1 x,float1 y,float1 z,float1 t,vec3 &E2,vec3 &B2)	//testhez
{
	buildTensorF_base(x,y,z,t);
	
	E2.x=tensorF.m[0][1]; //OK
	E2.y=tensorF.m[0][2];
	E2.z=tensorF.m[0][3];

	B2.x=tensorF.m[2][3];// inverse?  SI!   CHEAT!      (in lecture on Feynman -!)  +OK!!
	B2.y=tensorF.m[3][1];
	B2.z=tensorF.m[1][2];
}	
void buildBoostMatrixV()// tensorV
{
	vec3 n=normalize(move_dir);//move_dir_tr UA
	float1 y=gamma2, yi=gamma2-1.0;

	tensorV.m[0][0]=y;
	tensorV.m[0][1]=-y*beta*n.x;
	tensorV.m[0][2]=-y*beta*n.y;
	tensorV.m[0][3]=-y*beta*n.z;

	tensorV.m[1][0]=-y*beta*n.x;
	tensorV.m[1][1]=yi*n.x*n.x+1.0;
	tensorV.m[1][2]=yi*n.x*n.y;
	tensorV.m[1][3]=yi*n.x*n.z;

	tensorV.m[2][0]=-y*beta*n.y;
	tensorV.m[2][1]=yi*n.y*n.x;
	tensorV.m[2][2]=yi*n.y*n.y+1.0;
	tensorV.m[2][3]=yi*n.y*n.z;

	tensorV.m[3][0]=-y*beta*n.z;
	tensorV.m[3][1]=yi*n.z*n.x;
	tensorV.m[3][2]=yi*n.z*n.y;
	tensorV.m[3][3]=yi*n.z*n.z+1.0;

}			
//press 4
void fncdAtransfdA(float1 x,float1 y,float1 z,float1 t,vec3 &E2,vec3 &B2)	//origin ok   transform dA 
{
float1 t2=0.0;////                                   length contraction
vec3 P(x,y,z);vec3 P2=lorentz_transformationA(t,P,t2,LT1);x=P2.x;y=P2.y;z=P2.z; t=t2;  //field position!  2?3
//ugyan ugy mint a transformEMfield , trafo legelol 

	buildTensorF_base(x,y,z,t);
	
	buildBoostMatrixV();
qmatInverse(&tensorV,&tensorV); //test OK

	tensorVVFProduct();//!!!!!!!!!!!!!!!!!11
	
	E2.x=tensorF2.m[0][1];// A
	E2.y=tensorF2.m[0][2];
	E2.z=tensorF2.m[0][3];

	B2.x=tensorF2.m[2][3];//OK
	B2.y=tensorF2.m[3][1];
	B2.z=tensorF2.m[1][2];
}	
//press 5
void fncdAtransfA(float1 x,float1 y,float1 z,float1 t,vec3 &E2,vec3 &B2)	//test,    transform A
{
//float1 t2=0.0;//trafo itt nem jo, melyebben kell fncA4_tr
	buildTensorF_tr(x,y,z,t);//    ->dA
	
	E2.x=tensorF.m[0][1];// dA
	E2.y=tensorF.m[0][2];
	E2.z=tensorF.m[0][3];

	B2.x=tensorF.m[2][3];//OK
	B2.y=tensorF.m[3][1];
	B2.z=tensorF.m[1][2];
}	
		
vec3 lorentz_transformation4(float1 t1,vec3 P1,float1 &t2,int inv)
{						
	buildBoostMatrixV();
	if(inv&16) qmatInverse(&tensorV,&tensorV); //       KELL!!
/*	qmatTranspose(&tensorV,&tensorV);*/   // semmi NA

	if(inv&4) t1=-t1;

	float1 A42[4],A4[4]={t1,P1.x,P1.y,P1.z};

	for(int a=0;a<4;a++)//16 step 
	{
		A42[a]=0.0;
		for(int i=0;i<4;i++)
			A42[a] += (tensorV.m[a][i]*A4[i]);
//			A42[a] += (tensorV.m[i][a]*A4[i]);// pff    NA
	}
	vec3 P2;
	t2=A42[0];
	P2.x=A42[1];
	P2.y=A42[2];
	P2.z=A42[3];

	if(inv&4) t2=-t2;

	return P2;
}	
	
//--------------------------------------------------------------------------------------
//press 3
void transformEMfield(vec3 E,vec3 B,vec3 &E2,vec3 &B2 ,float1 x,float1 y,float1 z,float1 t)
{
	buildBoostMatrixV();

#if 0
//vector potential -bol          nyilvan jo
	buildTensorF_base(x,y,z,t);

#else
//E B -mezobol            ok       UA
	tensorF.m[0][0]=0.0;//++++            gaussian unit
	tensorF.m[0][1]=E.x/c;//c=1 soo        
	tensorF.m[0][2]=E.y/c;
	tensorF.m[0][3]=E.z/c;

	tensorF.m[1][0]=-E.x/c;//-++-
	tensorF.m[1][1]=0.0;
	tensorF.m[1][2]=B.z;              // BZ [1][2]  +!
	tensorF.m[1][3]=-B.y;

	tensorF.m[2][0]=-E.y/c;//--++
	tensorF.m[2][1]=-B.z;
	tensorF.m[2][2]=0.0;
	tensorF.m[2][3]=B.x;              // BX [2][3]

	tensorF.m[3][0]=-E.z/c;//-+-+
	tensorF.m[3][1]=B.y;              // BY [3][1]
	tensorF.m[3][2]=-B.x;
	tensorF.m[3][3]=0.0;

//NA	for(int a=0;a<4;a++)	for(int b=0;b<4;b++) if(a==0 || b==0) tensorF.m[a][b]*=-1.0;//upper lower indices
#endif

//KELL
qmatInverse(&tensorV,&tensorV);//test ?  NA //qmatTranspose(&tensorV,&tensorV);  OK         ugyan az csak nincs - benne lol


#if 0
//#ifdef _EDUMP
//	printf("E.z,dAw.z,dAz.w   %e %e %e  ,  %e %e\n",E.z,dAw.z,dAz.w,q1,q2);
	for(int a=0;a<4;a++)
	for(int b=0;b<4;b++) printf("%0.5f %0.5f %e\n",tensorF.m[a][b], tensordA.m[a][b], tensordA.m[a][b]/tensorF.m[a][b]);
	printf("\n");
#endif
	
	//derivalas utan transformalodik, szoval ben az A vektorpotential transformalodik,   csalas az egesz

		

	tensorVVFProduct();//!!!!!!!!!!!!!!!!!11
#if 0
	for(int a=0;a<4;a++)
	{
		for(int b=0;b<4;b++)  printf("%e ",tensorF2.m[a][b]);
		printf("\n");
	}
	printf("\n");
#endif	
	E2.x=tensorF2.m[0][1];
	E2.y=tensorF2.m[0][2];
	E2.z=tensorF2.m[0][3];

	B2.x=tensorF2.m[2][3];//OK +
	B2.y=tensorF2.m[3][1];
	B2.z=tensorF2.m[1][2];
}

//--------------------------------------------------------------------------------------

#if 0

F(ab) = dA(b)/dx(a) - dA(a)/dx(b) = 
= dAa.b - dAb.a
F[0][0]=dAt.w - dAt.w
F[0][1]=dAt.x - dAx.w
F[0][2]=dAt.y - dAy.w
F[0][3]=dAt.z - dAz.w

F[1][0]=dAx.w - dAt.x
F[1][1]=dAx.x - dAx.x
F[1][2]=dAx.y - dAy.x
F[1][3]=dAx.z - dAz.x

F[2][0]=dAy.w - dAt.y
F[2][1]=dAy.x - dAx.y
F[2][2]=dAy.y - dAy.y
F[2][3]=dAy.z - dAz.y

F[3][0]=dAz.w - dAt.z
F[3][1]=dAz.x - dAx.z
F[3][2]=dAz.y - dAy.z
F[3][3]=dAz.z - dAz.z





//Lorentz transformation - Wikipedia.html
const char *tensorF[]
{
/*
"0", "-Ex/c","-Ey/c","-Ez/c",// SI
"Ex/c",  "0", "-Bz",    "By",
"Ey/c", "Bz",   "0",   "-Bx",
"Ez/c","-By",  "Bx",     "0"
*/

  "0", "Ex", "Ey", "Ez",// gaussian unit
"-Ex",  "0", "Bz","-By",
"-Ey","-Bz",  "0", "Bx",
"-Ez", "By","-Bx",   "0"
};
const char *tensorV[]
{
  "y","-yB","0","0",  // y=gamma  B=beta=v/c
"-yB",  "y","0","0",
  "0",  "0","1","0",
  "0",  "0","0","1"
};
// X' = B(v)X

boost matrix is B(v)
    y,       -yBnx,       -yBny,       -yBnz,
-yBnx, 1+(y-1)nxnx,   (y-1)nxny,   (y-1)nxnz,
-yBny,   (y-1)nynx, 1+(y-1)nyny,   (y-1)nynz, 
-yBnz,   (y-1)nznx,   (y-1)nzny, 1+(y-1)nznz,
#endif


//--------------------------------------------------------------------------------------

void drawscene()
{
#define NN 1024
	float1 scale1=7.0*200.0*3;
	float1 scale2=scale1;// *0.5;

int h=32*3;
int s=32;//32 16

delta_e=1.0;
needarrowhead=01;//0

	gamma2=1.0/sqrt(1.0-v*v/(c*c));
	beta=v/c;

	
	float1 t2=0.0;
	vec3 v1,v2;
	
//move_dir_tr	
	v2=move_dir;//NA 
	move_dir_tr=normalize(lorentz_transformation2(0.0,v2,t2,0)/c);// NA  egyaltalan valtozik valami?

//Tdir_ty
	t2=0.0;
	v2=Tdir*200.0;
//vec3 testvec=normalize(v2-v1);	//!!!!!!!!
	v2=lorentz_transformation2(0.0,v2,t2,1)/c;// 1 OK ( no length contr) !!!!!!!!!!!!!!  1
	Tdir_tr=normalize(v2);//??????? ok?

//draw
	v2=Tdir*200.0;
	v2=lorentz_transformation2(0.0,v2,t2,0)/c;//0 OK
	arrow3d(v1,v2,0xffff00);
	v2=move_dir*200.0;
	arrow3d(v1,v2,0xff00ff);

//testvec
vec3 testvec=normalize(Tdir);// OK
	testvec=lorentz_transformation2(0.0,Tdir,t2,0)/c;// 0 OK
	arrow3d(v1,testvec*200.0,0x0000ff);
	
	vec3 q1=move_dir*v*gamma2;
	arrow3d(testvec*200.0,(testvec-q1)*200.0,0x00ff00);
	testvec -= move_dir*v*gamma2;// a lenyeg !!!!!!!!!!!!!!!!!!!!!!!!!!!
	arrow3d(v1,testvec*200.0,0xff0000);
	testvec=normalize(testvec);
	
	
	
	float1 maxval=-1e6;
	float1 minval=1e6;



if(viewmode==1||viewmode==4||viewmode==5)
{
	for(int z2=1;z2<NN;z2+=s)
	for(int y2=1;y2<NN;y2+=s)
	for(int x2=1;x2<NN;x2+=s)
//	if((z>=NN/2 -h)&&(z<=NN/2 +h) ||  (y==NN/2 ) ||	   (x==NN/2 ) )
	if((z2==NN/2+1 ) || 	(y2==NN/2+1) ||   (x2==NN/2+1) )
#ifdef _EandB
	if(x2>=NN/2)if(y2>=NN/2)if(z2>=NN/2) //OFF
#endif
	{
		float1 x=x2-512,y=y2-512,z=z2-512,t=0.0+time2;
		vec3 P1=vec3(x,y,z);

		if(espacetime) {y=0.0;t=y2-512;}

		int color=0x0000ff;
		if(viewmode==4) color=0x224477;
		if(viewmode==5) color=0x4400ff;

		vec3 E1,B1;	
		if(viewmode==1) E1=fncEgrad_tr(x,y,z,t);//view electric field    ???//rossz
		if(viewmode==4) fncdAtransfdA(x,y,z,t,E1,B1);// OK
		if(viewmode==5) fncdAtransfA(x,y,z,t,E1,B1);// OK

		
		float1 val=fabs(dot(normalize(E1),testvec));//merolegesseg teszt
		if(val>maxval) maxval=val;
		if(val<minval) minval=val;

		if(espacetime) {E1.y=0.0;}
		
		arrow3d(P1,P1+E1*scale1,color);//yellow blue
	}
	printf("E: %0.3f %0.3f \n",acos(minval)/iradian,acos(maxval)/iradian);
}


if(viewmode==1||viewmode==4||viewmode==5)//3
{
	 maxval=-1e6;
	 minval=1e6;

	for(int z2=1;z2<NN-1;z2+=s)
	for(int y2=1;y2<NN-1;y2+=s)
	for(int x2=1;x2<NN-1;x2+=s)
//	if((z>=NN/2 -h)&&(z<=NN/2 +h) || 	(y==NN/2+1) ||   (x==NN/2+1) )
	if((z2==NN/2+1 ) || 	(y2==NN/2+1) ||   (x2==NN/2+1) )
#ifdef _EandB
	if(x2>=NN/2)if(y2>=NN/2)if(z2>=NN/2) // OFF
#endif
	{
		float1 x=x2-512,y=y2-512,z=z2-512,t=0.0+time2;
		vec3 P1=vec3(x,y,z);

		if(espacetime) {y=0.0;t=y2-512;}

		int color=0xff0000;
		if(viewmode==4) color=0x774422;
		if(viewmode==5) color=0xff0044;

		vec3 E1,curl;	
		if(viewmode==1) curl=curlA_tr(x,y,z,t);// 
		if(viewmode==4) fncdAtransfdA(x,y,z,t,E1,curl);//OK
		if(viewmode==5) fncdAtransfA(x,y,z,t,E1,curl);;

		float1 val=fabs(dot(normalize(curl),testvec));//merolegesseg teszt
		if(val>maxval) maxval=val;
		if(val<minval) minval=val;

		if(espacetime) {curl.y=0.0;}
		
//		arrow3d(P1,P1+curl*scale1,getcolor(length(curl),0.9));//rainbow
		arrow3d(P1,P1+curl*scale1,color);//red       curl of A  after Lorentz transformation of A4 
	}

	printf("B: %0.3f %0.3f \n",acos(minval)/iradian,acos(maxval)/iradian);
	}	
	

if(viewmode==2)
{
	 maxval=-1e6;
	 minval=1e6;

	for(int z2=1;z2<NN;z2+=s)
	for(int y2=1;y2<NN;y2+=s)
	for(int x2=1;x2<NN;x2+=s)
//	if((z>=NN/2 -h)&&(z<=NN/2 +h) ||  (y==NN/2 ) ||	   (x==NN/2 ) )
	if((z2==NN/2+1 ) || 	(y2==NN/2+1) ||   (x2==NN/2+1 ) )
#ifdef _EandB
	if(x2>=NN/2)if(y2>=NN/2)if(z2>=NN/2)//OFF
#endif
	{
		float1 x=x2-512,y=y2-512,z=z2-512,t=0.0+time2;
		vec3 P1=vec3(x,y,z);

		if(espacetime) {y=0.0;t=y2-512;}
			
		vec3 A1=fncA_tr(x,y,z,t);//view vector potential
		if(espacetime)
		{
			vec4 A4=fncA4_tr(x,y,z,t); 
			A1.x=A4.x; 
			A1.y=A4.w; //  y=time
			A1.z=A4.z;
		}

		float1 val=fabs(dot(normalize(A1),testvec));//merolegesseg teszt  90
		if(val>maxval) maxval=val;
		if(val<minval) minval=val;

#if 0
		{//mindig null, kirva anyjat aki kitalalta
		vec4 A=fncA4_tr(x,y,z,t); 
		vec4 A2=A; A2.w=-A2.w;  //igy jo csak nyilvan   ab-ba=0
		//  A2.x=-A2.x;A2.y=-A2.y;A2.z=-A2.z; // igy is 0
//		printf("%e %e %e %e \n",A.x,A.y,A.z,A.w);//OK

		tensorF2.m[0][0]=A.w*A2.w - A2.w*A.w;//dump wedge
		tensorF2.m[0][1]=A.w*A2.x - A2.x*A.w; // AA2 AA2  jo de az faszsag
		tensorF2.m[0][2]=A.w*A2.y - A2.y*A.w;
		tensorF2.m[0][3]=A.w*A2.z - A2.z*A.w;
	
		tensorF2.m[1][0]=A.x*A2.w - A2.w*A.x;
		tensorF2.m[1][1]=A.x*A2.x - A2.x*A.x;
		tensorF2.m[1][2]=A.x*A2.y - A2.y*A.x;
		tensorF2.m[1][3]=A.x*A2.z - A2.z*A.x;
	
		tensorF2.m[2][0]=A.y*A2.w - A2.w*A.y;
		tensorF2.m[2][1]=A.y*A2.x - A2.x*A.y;
		tensorF2.m[2][2]=A.y*A2.y - A2.y*A.y;
		tensorF2.m[2][3]=A.y*A2.z - A2.z*A.y;
	
		tensorF2.m[3][0]=A.z*A2.w - A2.w*A.z;
		tensorF2.m[3][1]=A.z*A2.x - A2.x*A.z;
		tensorF2.m[3][2]=A.z*A2.y - A2.y*A.z;
		tensorF2.m[3][3]=A.z*A2.z - A2.z*A.z;

		float q=0.0;		
		for(int i=0;i<4;i++)//13 21 32
		for(int j=0;j<4;j++)
			if(!(i==1 && j==3))
			if(!(i==2 && j==1))
			if(!(i==3 && j==2))
				 q += tensorF2.m[i][j];
		
		if(q!=0.0)
		{
			printf("%.3f %.3f %.3f %.3f \n",tensorF2.m[0][0],tensorF2.m[0][1],tensorF2.m[0][2],tensorF2.m[0][3]);
			printf("%.3f %.3f %.3f %.3f \n",tensorF2.m[1][0],tensorF2.m[1][1],tensorF2.m[1][2],tensorF2.m[1][3]);
			printf("%.3f %.3f %.3f %.3f \n",tensorF2.m[2][0],tensorF2.m[2][1],tensorF2.m[2][2],tensorF2.m[2][3]);
			printf("%.3f %.3f %.3f %.3f \n",tensorF2.m[3][0],tensorF2.m[3][1],tensorF2.m[3][2],tensorF2.m[3][3]);
			printf("\n");	
		}
		}
#endif


//printf("%e \n",length(A1));		
		arrow3d(P1,P1+A1*scale1*1e-2,0xffff00);//yellow
	}
	printf("A: %0.3f %0.3f \n",acos(minval)/iradian,acos(maxval)/iradian);
	}	





if(viewmode==3)
{
	 maxval=-1e6;
	 minval=1e6;

	for(int z2=1;z2<NN;z2+=s)
	for(int y2=1;y2<NN;y2+=s)
	for(int x2=1;x2<NN;x2+=s)
	if((x2==NN/2+1) || 	(y2==NN/2+1) ||   (z2==NN/2+1) )
#ifdef _EandB
	if(x2>=NN/2)if(y2>=NN/2)if(z2>=NN/2)
#endif
	{
		float1 x=x2-512,y=y2-512,z=z2-512,t=0.0+time2;
		vec3 P1=vec3(x,y,z);

		if(espacetime) {y=0.0;t=y2-512;}

//E=-grad fi - dA / dt 	          A/dt==0 NOW
//                                                   length contraction
vec3 P(x,y,z);vec3 P2=lorentz_transformationA(t,P,t2,LT1);x=P2.x;y=P2.y;z=P2.z; t=t2;  //field position!  2?3
//a masik IR valtott, szoval valahogy kovetni kellene,  OK KELL

		vec3 E=fncEgrad_base(x,y,z,t);
		vec3 B=curlA_base(x,y,z,t);
//		fncdA(x,y,z,t,E,B); //UA      OK!

//		vec3 B3,E3;		fncdA(x,y,z,E3,B3);// UA  ok
#if 0
/*		printf("%e %e %e \n",E.x,E.y,E.z);		printf("%e %e %e \n",E3.x,E3.y,E3.z);
		printf("%e %e %e \n",B.x,B.y,B.z);		printf("%e %e %e \n\n",B3.x,B3.y,B3.z);*/
		
		printf("%e %e %e \n",E.x/E3.x,E.y/E3.y,E.z/E3.z);  //OK
		printf("%e %e %e \n\n",B.x/B3.x,B.y/B3.y,B.z/B3.z);
#endif		
		
		vec3 B2,E2;
#if 0
//ha minden t- akkor csak ez jo, cheatx!
		transformEMfield(E,B,E2,B2,x,y,z,t);//free direction
#else
E=-E;//KELL
//csak X-be !!!!!!!!
		B2.x=B.x;
		B2.y=(B.y+E.z*beta/c)*gamma2;//      X iranyu!!!!!!!!!!!!!!!!!!!!!!!!1
		B2.z=(B.z-E.y*beta/c)*gamma2;
		E2.x=E.x;
		E2.y=(E.y-B.z*beta)*gamma2;//B==0 szoval csak osszemegy!
		E2.z=(E.z+B.y*beta)*gamma2;//      *c????????
E2=-E2;		
#endif
		if(espacetime) {E2.y=0.0;B2.y=0.0;}
		
		arrow3d(P1,P1+B2*scale2,0xff4400);//    ordinary E/B Lorentz transformation
		arrow3d(P1,P1+E2*scale2,0x0044ff);//
		
		float1 val=fabs(dot(normalize(B2),testvec));//merolegesseg teszt  90
		if(val>maxval) maxval=val;
		if(val<minval) minval=val;
		val=fabs(dot(normalize(E2),testvec));//merolegesseg teszt  75
		if(val>maxval) maxval=val;
		if(val<minval) minval=val;
	}
	
	printf("EB: %0.3f %0.3f \n",acos(minval)/iradian,acos(maxval)/iradian);
}		


#ifdef _TESTCHARGE
//test charge
{//EONLY
	float1 t2=0.0;
//vec3 P1=vec3(512,250,250);
vec3 P1=vec3(0,250,250);
//P1=lorentz_transformation2(0.0,P1,t2,16);  //NA curlA_tr!

vec3 P2=P1,P1b=P1,P2b=P1;
vec3 V1=vec3(0.4,-0.4,-2.4);// - vec3(v,0.0,0.0); // c=1!!!!!!!
//v=-v;
V1.x=(V1.x+v)/(1.0+V1.x*v/(c*c));//ttps://en.wikipedia.org/wiki/Velocity-addition_formula
//v=-v;

for(int i=0;i<5000;i++)
{
	vec3 B,E;
	float1 t=i,	x=P1.x,	y=P1.y,	z=P1.z;
//t=0.0;	
	
	B=curlA_tr(x,y,z,t);
	E=fncEgrad_tr(x,y,z,t);
	P1+=V1;
//	V1 += E*2.0 + cross(V1,B)*20.0;
	V1 -= (E + cross(V1,B))*0.2;//0.2
//	V1 += cross(V1,B)*20.2;
	
//	P1b=lorentz_transformation2(t,P1,t2,2+16);  
	if(v==0.0)
	{
		v=-v_def;//-?
		P1b=lorentz_transformation2(t,P1,t2,0);  
		v=0.0;
	}
	else//0.8
	{
		v=v_def;
		P1b=lorentz_transformation2(t,P1,t2,0);  
		v=v_def;
	}
	
	arrow3d(P1,P2,0xffffff);
	arrow3d(P1b,P2b,0x00ffff);
	P2=P1;
	P2b=P1b;
}
}
#endif



//a mezo melle elmaszik
//#ifndef _EONLY
#ifdef _EandB
	time2+=1.01;
#endif

//E=-grad fi - dA / dt 	
//curl E = d curl A/dt = -dB/dt

}


vec4 fncA5(float1 t,float1 x,float1 y,float1 z)
{
	vec4 v4;
#if 0	
	v4.x=t*1.0; // csak x valtozik idoben
	v4.y=-z*2.0; // terben valtozik
	v4.z=x*2.0;
//	v4.z=z*2.0;//kiesik
#else
	float1 phase=z*k - t*k; // DXX 0.01   +time!!!! OK
	v4.x=sin(phase)*100.0;
	v4.y=0.0;
	v4.z= cos(phase)*100.0;  //  -+ OK !!! ha +time!   szamolashoz igy jo
	v4.w= cos(phase)*100.0;  // de  phase=-time soo ++, a 4FP elore mutat
#endif	
	return v4;
}
void testtensor()
{
	float1 x=0.152,y=0.624,z=0.2345,t=0.1234;
	vec4 dAx,dAy,dAz,dAw;
	float1 dt=DXX, dx=DXX, dy=DXX, dz=DXX;// -dt ? NO!
	x=frnd(1000.0);
	y=frnd(1000.0);
	z=frnd(1000.0);
	t=frnd(1000.0);
	
//four gradient   -+++
	dAw=(fncA5(t-dt,x   ,y   ,z   ) - fncA5(t+dt,x   ,y   ,z   ))/(dt*2.0);
	dAx=(fncA5(t   ,x+dx,y   ,z   ) - fncA5(t   ,x-dx,y   ,z   ))/(dx*2.0);
	dAy=(fncA5(t   ,x   ,y+dy,z   ) - fncA5(t   ,x   ,y-dy,z   ))/(dy*2.0);
	dAz=(fncA5(t   ,x   ,y   ,z+dz) - fncA5(t   ,x   ,y   ,z-dz))/(dz*2.0);

//exterior derivative  Aij - Aji
	tensorF.m[0][0]= (dAw.w - dAw.w);// w == time
	tensorF.m[0][1]= (dAw.x - dAx.w);//Ex   idoderival - scalar derivaltja
	tensorF.m[0][2]= (dAw.y - dAy.w);//Ey
	tensorF.m[0][3]= (dAw.z - dAz.w);//Ez

	tensorF.m[1][0]= (dAx.w - dAw.x);
	tensorF.m[1][1]=-(dAx.x - dAx.x);
	tensorF.m[1][2]=-(dAx.y - dAy.x);//Bz
	tensorF.m[1][3]=-(dAx.z - dAz.x);

	tensorF.m[2][0]= (dAy.w - dAw.y);
	tensorF.m[2][1]=-(dAy.x - dAx.y);
	tensorF.m[2][2]=-(dAy.y - dAy.y);
	tensorF.m[2][3]=-(dAy.z - dAz.y);//Bx

	tensorF.m[3][0]= (dAz.w - dAw.z);
	tensorF.m[3][1]=-(dAz.x - dAx.z);//By  elvileg ([0][0] [1][2] [2][3] [3][1] egy 4d CURL , es 4 van belole)
	tensorF.m[3][2]=-(dAz.y - dAy.z);
	tensorF.m[3][3]=-(dAz.z - dAz.z);

//	if(fabs(tensorF.m[0][3])>1e-10 || fabs(tensorF.m[1][2])>1e-10)
if(1)
	{
		printf("E %.3f %.3f %e \n",tensorF.m[0][1],tensorF.m[0][2],tensorF.m[0][3]);
		printf("B %.3f %.3f %e \n",tensorF.m[2][3],tensorF.m[3][1],tensorF.m[1][2]);
		printf("\n");
	}
	else printf("OK\n");
#if 0	
	// z fuggese z tol eltunik,  kiszedte!
	printf("%.3f %.3f %.3f %.3f \n",dAw.w,dAw.x,dAw.y,dAw.z);
	printf("%.3f %.3f %.3f %.3f \n",dAx.w,dAx.x,dAx.y,dAx.z);
	printf("%.3f %.3f %.3f %.3f \n",dAy.w,dAy.x,dAy.y,dAy.z);
	printf("%.3f %.3f %.3f %.3f \n",dAz.w,dAz.x,dAz.y,dAz.z);
	printf("\n");
	printf("%.3f %.3f %.3f %.3f \n",tensorF.m[0][0],tensorF.m[0][1],tensorF.m[0][2],tensorF.m[0][3]);
	printf("%.3f %.3f %.3f %.3f \n",tensorF.m[1][0],tensorF.m[1][1],tensorF.m[1][2],tensorF.m[1][3]);
	printf("%.3f %.3f %.3f %.3f \n",tensorF.m[2][0],tensorF.m[2][1],tensorF.m[2][2],tensorF.m[2][3]);
	printf("%.3f %.3f %.3f %.3f \n",tensorF.m[3][0],tensorF.m[3][1],tensorF.m[3][2],tensorF.m[3][3]);
	printf("\n");

	for(int a=0;a<4;a++)// Fab
	for(int b=0;b<4;b++) 
	{
		tensorF3.m[a][b]=tensorF.m[a][b];
		if(a==0 || b==0) tensorF3.m[a][b]=-tensorF.m[a][b];//  -E
	}
	
	for(int a=0;a<4;a++)// F^ab * Fab   nem transpose ,     Fab=F^ab de timecomponent elojelet valt (E)
	for(int b=0;b<4;b++)
	{
		tensorF2.m[a][b]=tensorF.m[a][b]*tensorF3.m[a][b];//negyzetre emeli, es - minden E (time components)
	}
	printf("%.3f %.3f %.3f %.3f \n",tensorF2.m[0][0],tensorF2.m[0][1],tensorF2.m[0][2],tensorF2.m[0][3]);
	printf("%.3f %.3f %.3f %.3f \n",tensorF2.m[1][0],tensorF2.m[1][1],tensorF2.m[1][2],tensorF2.m[1][3]);
	printf("%.3f %.3f %.3f %.3f \n",tensorF2.m[2][0],tensorF2.m[2][1],tensorF2.m[2][2],tensorF2.m[2][3]);
	printf("%.3f %.3f %.3f %.3f \n",tensorF2.m[3][0],tensorF2.m[3][1],tensorF2.m[3][2],tensorF2.m[3][3]);
	printf("\n");	
	
	
	
//exterior derivative , 
	vec4 A=fncA5(t,x ,y ,z );
	
	tensorF2.m[0][0]=A.w*A.w - A.w*A.w;
	tensorF2.m[0][1]=A.w*A.x - A.x*A.w;
	tensorF2.m[0][2]=A.w*A.y - A.y*A.w;
	tensorF2.m[0][3]=A.w*A.z - A.z*A.w;

	tensorF2.m[1][0]=A.x*A.w - A.w*A.x;
	tensorF2.m[1][1]=A.x*A.x - A.x*A.x;
	tensorF2.m[1][2]=A.x*A.y - A.y*A.x;
	tensorF2.m[1][3]=A.x*A.z - A.z*A.x;

	tensorF2.m[2][0]=A.y*A.w - A.w*A.y;
	tensorF2.m[2][1]=A.y*A.x - A.x*A.y;
	tensorF2.m[2][2]=A.y*A.y - A.y*A.y;
	tensorF2.m[2][3]=A.y*A.z - A.z*A.y;

	tensorF2.m[3][0]=A.z*A.w - A.w*A.z;
	tensorF2.m[3][1]=A.z*A.x - A.x*A.z;
	tensorF2.m[3][2]=A.z*A.y - A.y*A.z;
	tensorF2.m[3][3]=A.z*A.z - A.z*A.z;
	
	printf("%.3f %.3f %.3f %.3f \n",tensorF2.m[0][0],tensorF2.m[0][1],tensorF2.m[0][2],tensorF2.m[0][3]);
	printf("%.3f %.3f %.3f %.3f \n",tensorF2.m[1][0],tensorF2.m[1][1],tensorF2.m[1][2],tensorF2.m[1][3]);
	printf("%.3f %.3f %.3f %.3f \n",tensorF2.m[2][0],tensorF2.m[2][1],tensorF2.m[2][2],tensorF2.m[2][3]);
	printf("%.3f %.3f %.3f %.3f \n",tensorF2.m[3][0],tensorF2.m[3][1],tensorF2.m[3][2],tensorF2.m[3][3]);
	printf("\n");	
	
#endif	
	
	return;
}



int main()
{

//for(int i=0;i<100;i++) testtensor();return 0;

	init_system();
	
	int first=1;
	float1 alpha=30,beta=30;


	while(1)
	{
		if(keyboard() || first||key)
		{
			first=0;

			qclear();
//camera		
//		eye=vec3(605,500,300);
			eye=sphere(alpha,beta)*1000.0;
			look=vec3(0,110,0);
			eye+=look;
			up=vec3(0,1,0);

			alpha+=dmx*0.6;
			beta+=dmy*0.6;

		
			if(key=='q') v=0.0*c;
			if(key=='w') v=v_def*c;
//			if(key=='t') etradx^=1;
			if(key==' ') espacetime^=1;

			

			calc_camaxis();
			drawscene();

	
			qflush();//		qsleep(10);
			key=0;
		}
	}
//	getchar();
return 0;
}

//vec3 P(x,y,z);vec3 P2=lorentz_transformation2(t,P,t2);x=P2.x;y=P2.y;z=P2.z; t=t2;  //field position!  2?3

//NA
#if 0
void transformDX(float1 &dx,float1 &dy,float1 &dz,float1 &dt)
{
	float1 t2=0.0;
	vec3 P=vec3(dx,dy,dz);
	vec3 P2=lorentz_transformation2(dt,P,t2); //field delta position!    2        2i & 3  UA!
	
	dx=P2.x;
	dy=P2.y;
	dz=P2.z; 
	dt=t2; 
}
//	float1 q=DXX-DXX/gamma2;
//	float1 dx=DXX- q*(fabs( move_dir.x));//  ? sqr sqrt                        EZ ROSSZ!!!!!!!!!!!!!!!
#endif

//???????wedge product       + [Ai,Aj]= AiAj - AjAi         ([C,D]= CD - DC )
