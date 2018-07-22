
#define __USE_X11


#include "eng6.cpp"


vec3 viewmatx,viewmaty,viewmatz;

void lineproj(vec3 p1,vec3 p2,int col)
{
    int x1[2],y1[2];

    x1[0]=dot(p1,viewmatx) +1200/2; 
    y1[0]=dot(p1,viewmaty) +600/2;
//  z1[0]=dot(p1,viewmatz);

    x1[1]=dot(p2,viewmatx) +1200/2; 
    y1[1]=dot(p2,viewmaty) +600/2;
//  z1[1]=dot(p2,viewmatz);

	if(length(p1-p2)<500.0)
	    line(x1[0],y1[0],x1[1],y1[1],col);
}


vec3 fnc(float1 u,float1 v)//getsurfacepoint
{
	float1 rs=1.0;
	vec3 v1;
	float1 r=1.0/(sqrt(u*u+v*v)*0.001);
	v1.x=u;
	v1.y=2.0*sqrt(rs*(r-rs))*40.0;//Flamm*40
	v1.z=v;
	return v1;
}
vec3 getsurfacenormal(vec3 p0)
{
	float1 u=p0.x;
	float1 v=p0.z;
	
	vec3 p1=fnc(u,v);
	vec3 p2=fnc(u+0.1,v);
	vec3 p3=fnc(u,v+0.1);
	
	return normalize(cross(p2-p1,p3-p1));
}
int main()
{
	init_system();

	float1 a1=30.0*iradian,b1=20.0*iradian;
	viewmatz=vec3(sin(a1)*cos(b1),sin(b1),cos(a1)*cos(b1));
	viewmatx=normalize(cross(viewmatz,vec3(0,1,0)));
	viewmaty=normalize(cross(viewmatx,viewmatz));

//surface
	int h=20;
	for(int y=-600;y<600;y+=h)
	for(int x=-600;x<600;x+=h)
	{
	    vec3 d1=fnc(x  ,y);
	    vec3 d2=fnc(x+h,y);
	    vec3 d3=fnc(x  ,y+h);

	    lineproj(d1,d2,0x00aa00);
    	lineproj(d1,d3,0x00aa00);
	}
//moving point
//	for(int a=0;a<5;a++)
	int a=0;
	{
	float1 wavelength=0.0;
	if(a==0) wavelength=0.2;
	if(a==1) wavelength=2.0;
	if(a==2) wavelength=5.0;
	if(a==3) wavelength=10.0;
	if(a==4) wavelength=30.0;
	
	vec3 pos=fnc(500,80);
	vec3 pos2=pos;
	vec3 vel(-20,0,0);
	float1 alpha=150.0*iradian;
	if(a) vel=vec3(
		-20*sin(90.0*iradian-alpha/2.0),
		0,
		-20*cos(90.0*iradian-alpha/2.0));
	
	int o=0;
	float1 dt=0.01,s=0;
	for(int i=0;i<100000;i++)
	{
	    lineproj(pos,pos2,color(255,250,0));//a

		vec3 normal=getsurfacenormal(pos);
		pos+=vel*dt;
		vel-=normal*dot(normal,vel);
#if 0		
		if(a!=0)
		if(s>=wavelength) 
		{
			vec3 bitangent=normalize(cross(vel,normal));
/*			
			if(o==0) vel=bitangent*length(vel);
			else     vel=bitangent*(-length(vel));*/
			float1 alpha2=alpha;
			if(o==1) alpha2=-alpha2;
			vel=(bitangent*sin(alpha2) + normalize(vel)*cos(alpha2))*length(vel);
			o^=1;
			s=0;
		}
#endif
		s+=length(vel)*dt;
		
		pos2=pos;
	}
	}
	
	qflush();
	getchar();

	return 0;
}



