
#define __USE_X11

#include "eng6.cpp"


//#define SPHERE
#define TORUS
//#define SURFACE

//#define __VIEW3D


#define pi M_PI




int iter=30;//15


vec3 P(float1 u,float1 v,float1 R)
{
    vec3 v_;

#ifdef SPHERE
//sphere
    v_.x=R*cos(u)*sin(v);
    v_.y=R*sin(u)*sin(v);
    v_.z=R*cos(v);//  sphere
//    v_.z=R*cos(v)/1.3;// flat
#endif
#ifdef TORUS
    float1 r2=50.0;
    R+=r2*cos(v);//R=c   r2=a  (c+a*cos(v))=R
    v_.x=R*cos(u);
    v_.y=R*sin(u);//R=a
    v_.z=r2*sin(v);//r2=a
#endif
#ifdef SURFACE
    v_.x=(u-0.5)*100.0;
    v_.y=20.0*sin(u*4.0)*sin(v*3.0);
    v_.z=(v-0.5)*100.0;
#endif

    return v_;
}
vec3 get_surface_normal(vec3 pos,float1 R)
{
#ifdef SPHERE
    return normalize(pos);
#endif
#ifdef TORUS
    vec3 axis=normalize(vec3(pos.x,pos.y,0.0f))*R;
    
    return normalize(pos-axis);
#endif
#ifdef SURFACE
#endif
}
vec3 get_surface_normalUV(float1 u,float1 v,float1 R)
{
#ifdef SPHERE
#endif
#ifdef TORUS
    vec3 p1,p2;

    p1.x=R*cos(u);
    p1.y=R*sin(u);//R=a
    p1.z=0.0;

    float1 r2=50.0;
    R+=r2*cos(v);//R=c   r2=a  (c+a*cos(v))=R
    p2.x=R*cos(u);
    p2.y=R*sin(u);//R=a
    p2.z=r2*sin(v);//r2=a  
    
    return normalize(p2-p1);
#endif
#ifdef SURFACE
#endif
}


vec3 derivative_ofP(float1 u,float1 v,float1 du,float1 dv,float1 R)
{
    return (P(u+du,v+dv,R) - P(u,v,R)); // /dx !!!!!!! 
}


vec3 viewmatx,viewmaty,viewmatz;


void lineproj(vec3 pos,vec3 p2,int col)
{
    float1 x1[2],y1[2];

    x1[0]=dot(pos,viewmatx) +float1(1200/2); 
    y1[0]=dot(pos,viewmaty) +float1(600/2);
//  z1[0]=dot(pos,viewmatz);

    x1[1]=dot(p2,viewmatx) +float1(1200/2); 
    y1[1]=dot(p2,viewmaty) +float1(600/2);
//  z1[1]=dot(p2,viewmatz);

    line(x1[0] ,y1[0] ,x1[1] ,y1[1] ,col);
}


void line2d(vec3 vel,int col)
{
	line(300,300,300+vel.x ,300+vel.y ,col);
}


void grav()// Bernoulli grav
{
	float1 a,m=5.8e24,c=3e8,g=6.7e-11,r=6370e3,m2=10.0,F,rs,ro,A;
	
	rs=2.0*m*g/(c*c);

	r=rs*2.0;
	F=m*m*g/(r*r);
	
//	p= -F/A;    ro=m/V;     p= -v^2/2 *ro     F= -c^2/2 *ro  
	A=r*r;// A=1 m^2
	ro=(F/A)/(c*c/2.0);
	printf("%e \n",ro);
	
	m2=ro;
	r=2.0*rs;	m2=ro/(r*r*r);
	printf("%e \n",m2);
}


int main()
{
//grav();return 0;
init_system();

qclear();


float1 u_=23.0;
float1 v_=53.0;
// u_=63.0; v_=13.0;//test
//  u_=123.0; v_=73.0;//test   ok mind xD
float1 du_=-0.025;
float1 dv_=0.074;

int ii=1525000;
float1 u=u_*iradian;
float1 v=v_*iradian;
float1 du=du_*iradian;
float1 dv=dv_*iradian;

float1 z=0.0001*iradian;//0.0001
float1 d=0.0001*iradian;
float1 R=200;
float1 g;

float1 dt=0.05;//0.01  ii/=10;
float1 duu_du,duv_du,dvu_du,dvv_du,duu_dv,duv_dv,dvu_dv,dvv_dv;



{
float1 a1=30.0*iradian,b1=20.0*iradian;
viewmatz=vec3(sin(a1)*cos(b1),sin(b1),cos(a1)*cos(b1));
viewmatx=normalize(cross(viewmatz,vec3(0,1,0)));
viewmaty=normalize(cross(viewmatx,viewmatz));
}


//draw torus
for(int y=0;y<30;y++)
for(int x=0;x<30;x++)
{
    float1 u1=((float1)(x)/30.0)*M_PI*2;
    float1 vel=((float1)(y)/30.0)*M_PI*2;
    float1 u2=((float1)(x+1)/30.0)*M_PI*2;
    float1 v2=((float1)(y+1)/30.0)*M_PI*2;

    vec3 d1=P(u1,vel,R);
    vec3 d2=P(u1,v2,R);
    vec3 d3=P(u2,vel,R);

#ifdef __VIEW3D
    lineproj(d1,d2,0x550000);
    lineproj(d1,d3,0x550000);
#else
    line(300+(int)d1.x,300+(int)d1.y,300+(int)d2.x,300+(int)d2.y,0x550000);
    line(300+(int)d1.x,300+(int)d1.y,300+(int)d3.x,300+(int)d3.y,0x550000);

    line(800+(int)d1.x,300+(int)d1.z,800+(int)d2.x,300+(int)d2.z,0x550000);
    line(800+(int)d1.x,300+(int)d1.z,800+(int)d3.x,300+(int)d3.z,0x550000);
#endif
}
qflush();



//solution 1
//cut the normal component of velocity (simple parallel transport)
vec3 pos=P(u,v,R);
vec3 vel=P(u+du,v+dv,R) - P(u,v,R);


for(int i=0;i<ii;i++)
{
    vec3 normal=get_surface_normal(pos,R);

    pos=pos+vel*dt;
    vel=vel-normal*dot(normal,vel);

#ifdef __VIEW3D
    lineproj(pos,pos+vec3(1,1,1),0x00aaff);
#else
    pixel(300+(int)pos.x,300+(int)pos.y,0x00aaff);//light blue
    pixel(800+(int)pos.x,300+(int)pos.z,0x00aaff);
#endif
}
qflush();


//if(0)
{
float1 Cu_uu,Cu_uv,Cu_vu,Cu_vv,Cv_uu,Cv_uv,Cv_vu,Cv_vv;

//Parallel transport
 u=u_*iradian;
 v=v_*iradian;
 du=du_*iradian;
 dv=dv_*iradian;
pos=P(u,v,R);
vel=P(u+du,v+dv,R) - P(u,v,R);
    vec3 normal=get_surface_normal(pos,R);
    float1 curv;
    
    du=dv=0.0;
    
    vec3 Puu, Puv, Pvu, Pvv;
	vec3 Pu =derivative_ofP(u,v, z,0.0,R)/z;//first derivatives
	vec3 Pv =derivative_ofP(u,v, 0.0,z ,R)/z;
	vec3 Pu_=Pu/dot(Pu,Pu);//contravatiant of tangents
	vec3 Pv_=Pv/dot(Pv,Pv);
    
	for(int k=0;k<2;k++)
	for(int l=0;l<2;l++)
	{
		float1 dx=z;//0.0001*iradian; // z & d
		if(l==0)	{du=dx; dv=0.0;}
		if(l==1)	{dv=dx; du=0.0;}
		vec3 trans_vec=derivative_ofP(u,v,du,dv,R)/dx;
		
		float1 dalf=0.1*iradian; // less more precise
		float1 du2,dv2;
		
		int nn=1000;	
		for(int i=0;i<=nn;i++)
		{
			float1 t1=(float1)i/(float1)nn;
			if(k==0)	{du2=dalf*t1; dv2=0.0;}
			if(k==1)	{dv2=dalf*t1; du2=0.0;}
#if 0
			vec3 pos2=P(u+du2, v+dv2, R);//		vec3 vel=P(u+du,v+dv,R) - P(u,v,R); // OK
	    	vec3 normal2=get_surface_normal(pos2,R);
#else
	    	vec3 normal2=get_surface_normalUV(u+du2,v+dv2,R);// ua OK
		    trans_vec=trans_vec - normal2*dot(normal2,trans_vec);// essence of parallel transport ! = cutting the normal component of vector
#endif
//	NA	    trans_vec=normalize(trans_vec)*len;// scale   , NA :length is wrong!!
		}

		if(k==0)	{du2=dalf; dv2=0.0;}//just for sure
		if(k==1)	{dv2=dalf; du2=0.0;}
		vec3 local_tan=derivative_ofP(u+du2,v+dv2,du,dv,R)/dx;// at NEW POS!
	    vec3 delta_vec=(trans_vec - local_tan)/dalf;  // inverse second derivative (derivative of tangent)
   
	    du=dot(-delta_vec,Pu_); //  - ! because the parallel transport inverse of the second derivative
	    dv=dot(-delta_vec,Pv_);
	    
	    if(l==0) if(k==0) {Cu_uu=du;Cv_uu=dv;}
	    if(l==0) if(k==1) {Cu_uv=du;Cv_uv=dv;}
	    if(l==1) if(k==0) {Cu_vu=du;Cv_vu=dv;}
	    if(l==1) if(k==1) {Cu_vv=du;Cv_vv=dv;}
    }
    
	printf("%.3f   \n",curv);
//	printf("%.3f %.3f %.3f  \n",refvec.x,refvec.y,refvec.z);
	printf("%f %f %f %f \n",Cu_uu,Cu_uv,Cu_vu,Cu_vv);
	printf("%f %f %f %f \n",Cv_uu,Cv_uv,Cv_vu,Cv_vv);
	printf("\n");
}	
qflush();


	

//soulution 2
//simple function of second derivative
//http://www.win.tue.nl/~rvhassel/Onderwijs/Tensor-ConTeX-Bib/Examples-diff-geom/Torus-diff-geom/torus-together.pdf
 u=u_*iradian;
 v=v_*iradian;
 du=du_*iradian;
 dv=dv_*iradian;
pos=P(u,v,R);
vel=P(u+du,v+dv,R) - P(u,v,R);

for(int i=0;i<ii;i++)
{
	float1 Cu_uu,Cu_uv,Cu_vu,Cu_vv,Cv_uu,Cv_uv,Cv_vu,Cv_vv;
	float1 a=50.0,c=R;//a==r2

	Cu_uu=0.0;
	Cu_uv=(-a*sin(v))/(c+a*cos(v));
	Cu_vu=Cu_uv;
	Cu_vv=0.0;

	Cv_uu=(1.0/a)*(sin(v)*(c+a*cos(v)));
	Cv_uv=0.0;
	Cv_vu=0.0;
	Cv_vv=0.0;
	
	float1 ddu= -Cu_uu*du*du - Cu_uv*du*dv - Cu_vu*dv*du - Cu_vv*dv*dv;
	float1 ddv= -Cv_uu*du*du - Cv_uv*du*dv - Cv_vu*dv*du - Cv_vv*dv*dv;

	u += du*dt;
	v += dv*dt;
	du += ddu*dt;
	dv += ddv*dt;

	pos=P(u,v,R);

#ifdef __VIEW3D
    lineproj(pos,pos+vec3(1,1,1),0xff00ff);//purple
#else
    pixel(300+(int)pos.x,300+(int)pos.y,0xff00ff);
    pixel(800+(int)pos.x,300+(int)pos.z,0xff00ff);
#endif

	if(i==0)
	{
		printf("%f %f %f %f \n",Cu_uu,Cu_uv,Cu_vu,Cu_vv);
		printf("%f %f %f %f \n",Cv_uu,Cv_uv,Cv_vu,Cv_vv);
		printf("\n");
	}
}
qflush();



//solution 3-4
//makes derivative by vectors   (original Christoffel symbols)
for(int q=0;q<2;q++)
{
 u=u_*iradian;
 v=v_*iradian;
 du=du_*iradian;
 dv=dv_*iradian;

pos=P(u,v,R);
vel=P(u+du,v+dv,R) - P(u,v,R);


for(int i=0;i<ii;i++)
{
vec3 Pu = derivative_ofP(u,v, z,0.0,R)/z;//first derivatives
vec3 Pv = derivative_ofP(u,v, 0.0,z ,R)/z;
//printf("%e %e \n",length(Pu),length(Pv));

vec3 Pu_=Pu/dot(Pu,Pu);//contravatiant of tangents
vec3 Pv_=Pv/dot(Pv,Pv);


vec3 Pu_su = derivative_ofP(u+d,v, z,0.0,R)/z;// (s)hifted values
vec3 Pu_sv = derivative_ofP(u,v+d, z,0.0,R)/z;
vec3 Pv_su = derivative_ofP(u+d,v, 0.0,z,R)/z;
vec3 Pv_sv = derivative_ofP(u,v+d, 0.0,z,R)/z;


float1 guu = dot(Pu,Pu);//E
float1 guv = dot(Pu,Pv);//F   = float1 gvu = dot(Pv,Pu);
float1 gvv = dot(Pv,Pv);//G


vec3 Puu=(Pu_su-Pu)/d;//second derivatives
vec3 Puv=(Pu_sv-Pu)/d;
vec3 Pvu=(Pv_su-Pv)/d;
vec3 Pvv=(Pv_sv-Pv)/d;


float1 Cu_uu,Cu_uv,Cu_vu,Cu_vv,Cv_uu,Cv_uv,Cv_vu,Cv_vv;

if(q==0)//classic way
{
	duu_du = dot(Puu,Pu) + dot(Pu,Puu);//derivatives of components of metric tensor
	duu_dv = dot(Puv,Pu) + dot(Pu,Puv);
	duv_du = dot(Pvu,Pu) + dot(Pv,Puu);
	duv_dv = dot(Pvv,Pu) + dot(Pv,Puv);

	dvu_du = dot(Puu,Pv) + dot(Pu,Pvu);
	dvu_dv = dot(Puv,Pv) + dot(Pu,Pvv);
	dvv_du = dot(Pvu,Pv) + dot(Pv,Pvu);
	dvv_dv = dot(Pvv,Pv) + dot(Pv,Pvv);


	g=1.0/(2.0*guu);
	Cu_uu = (duu_du + duu_du - duu_du )*g;//Christoffel symbols (original form)
	Cu_uv = (duu_dv + duv_du - duv_du )*g;
	Cu_vu = (duv_du + duu_dv - dvu_du )*g;
	Cu_vv = (duv_dv + duv_dv - dvv_du )*g;

	g=1.0/(2.0*gvv);
	Cv_uu = (dvu_du + dvu_du - duu_dv )*g;
	Cv_uv = (dvu_dv + dvv_du - duv_dv )*g;
	Cv_vu = (dvv_du + dvu_dv - dvu_dv )*g;
	Cv_vv = (dvv_dv + dvv_dv - dvv_dv )*g;
}
else//simple way           myway!
{
	Cu_uu = dot(Puu,Pu_);//Christoffel symbols (simplified Gauss–Codazzi solution)
	Cu_uv = dot(Puv,Pu_);
	Cu_vu = dot(Pvu,Pu_);
	Cu_vv = dot(Pvv,Pu_);

	Cv_uu = dot(Puu,Pv_);
	Cv_uv = dot(Puv,Pv_);
	Cv_vu = dot(Pvu,Pv_);
	Cv_vv = dot(Pvv,Pv_);
}
	
	if(i==0)
	{
		printf("%f %f %f %f \n",Cu_uu,Cu_uv,Cu_vu,Cu_vv);
		printf("%f %f %f %f \n",Cv_uu,Cv_uv,Cv_vu,Cv_vv);
		printf("\n");
	}

float1 ddu= -Cu_uu*du*du - Cu_uv*du*dv - Cu_vu*dv*du - Cu_vv*dv*dv;
float1 ddv= -Cv_uu*du*du - Cv_uv*du*dv - Cv_vu*dv*du - Cv_vv*dv*dv;

u += du*dt;
v += dv*dt;
du += ddu*dt;
dv += ddv*dt;

pos=P(u,v,R);

int col=0x00ff00;//green
if(q) col=0xffff00;//yellow
#ifdef __VIEW3D
    lineproj(pos,pos+vec3(1,1,1),col);
#else
    pixel(300+(int)pos.x,300+(int)pos.y,col);
    pixel(800+(int)pos.x,300+(int)pos.z,col);
#endif
}
qflush();
}

qflush();
getchar();

return 0;
}



