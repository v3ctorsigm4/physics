
#define __USE_X11

#include "eng6.cpp"



float1 rs,dt,rs2,g,c,m,r1,v1,
    c2 = 1.0,
    scale2 = 1.0,
    scale;
vec3 pos, vel;
int n;
int dxx=0;


void pixel5(vec3 v1,int col)
{
    v1.x*=scale2;
    v1.y*=scale2;
    v1.z*=scale2;
    
//    pixel(500+v1.x+dxx,300+v1.y,col);
    XSetForeground(display5,gc,col);    
    XDrawPoint(display5, win, gc, 500+v1.x+dxx,300+v1.y);
}

void schwarzschild()
{
    float1 r,p,t,T, dr,dp, ddr,ddT,ddp,  dT;
    
    scale = c;//  c=1
    c2  = c /scale;
    rs2 = rs/scale;

    pos /= scale;
    vel /= scale;


    r  = sqrt(pos.x*pos.x + pos.y*pos.y);
    p  = atan2(pos.y,pos.x);
    t = 0.0;

    dr = (pos.x * vel.x + pos.y * vel.y) / r;
    dp = (pos.x * vel.y - pos.y * vel.x) / (r*r);
    dT=1.0;

    ddr=0;
    ddp=0;
    ddT=0;
    
    T = 0.0;
    


    for(int j = 0;j<n;j++)
    {
        if (r < rs2) return;
    
        float1 u=rs2/2.0;
        ddT  = -2.0*u/(r*(r-2.0*u)) * dT*dr;
        ddr  = -(u*(r-2.0*u)/(r*r*r))*dT*dT
                +(u/(r*(r-2.0*u)))*dr*dr
                +(r-2.0*u)*dp*dp;
        ddp  = -dr*dp*2.0/r;

    
        r+=dr*dt;     
        p+=dp*dt;     
        t+=dt;         

        dr+=ddr*dt;
        dp+=ddp*dt;
        dT+=ddT*dt;

        T+=dt/dT;    
        

        pos.x = r*cos(p)*scale;
        pos.y = r*sin(p)*scale;

        
        pixel5(pos,0x00ff00);
    }
}
void schwarzschild2()
{
    float1 r,p,t,T, dr,dp, ddr,ddT,ddp,  dT;
    
    scale = c;//  c=1
    c2  = c /scale;
    rs2 = rs/scale;

    pos /= scale;
    vel /= scale;


    r  = sqrt(pos.x*pos.x + pos.y*pos.y);
    p  = atan2(pos.y,pos.x);
    t = 0.0;

    dr = (pos.x * vel.x + pos.y * vel.y) / r;
    dp = (pos.x * vel.y - pos.y * vel.x) / (r*r);
    dT=1.0;

    ddr=0;
    ddp=0;
    ddT=0;
    
    T = 0.0;
    


    for(int j = 0;j<n;j++)
    {
        if (r < rs2) return;
        float1 u=rs2/2.0;

//ds^2=((r-2u)/r) dt^2 - (r/(r-2u)) dr^2 - (r^2) do^2 -(r^2 sin2 o) df^2

        vec3 Rt(sqrt((r-2.0*u)/r),0.0,0.0);//first derivative
        vec3 Rr(0.0,sqrt(r/(r-2.0*u)),0.0);
        vec3 Rp(0.0,0.0,sqrt(r*r));

		float1 gtt=dot(Rt,Rt);//metric tensor
		float1 grr=dot(Rr,Rr);
		float1 gpp=dot(Rp,Rp);
		
//		float1 grt=dot(Rr,Rt);//

/*		vec3 Rt_=Rt/gtt;//contravariant
		vec3 Rr_=Rr/grr;
		vec3 Rp_=Rp/gpp;*/


        float1 dx=1e-12;
        float1 r2=r+dx;
        vec3 Rt_dx(sqrt((r2-2.0*u)/r2),0.0,0.0);//neightb
        vec3 Rr_dx(0.0,sqrt(r2/(r2-2.0*u)),0.0);
        vec3 Rp_dx(0.0,0.0,sqrt(r2*r2));
        
        
		vec3 Rtr=(Rt_dx -Rt)/dx;//second derivative respect to r
		vec3 Rrr=(Rr_dx -Rr)/dx;
		vec3 Rpr=(Rp_dx -Rp)/dx;


        float1 Ttrt=dot(Rtr,Rt)/gtt;
        float1 Trtt=dot(Rtr,Rt)/grr;//no Rr'   dot(Rtt,Rr)/grr
        float1 Trrr=dot(Rrr,Rr)/grr;//change!  dot(Rrr,Rr)/grr
        float1 Trpp=dot(Rpr,Rp)/grr;// 		   dot(Rpp,Rr)/grr
        float1 Tprp=dot(Rpr,Rp)/gpp;

//        float1 Tttt=(dot((Rt_dx -Rt)/dx,Rt) +dot(Rt,(Rt_dx -Rt)/dx)) /(gtt*2.0);
        //float1 Trtt=(dot((Rt_dx -Rt)/dx,Rt) +dot(Rt,(Rt_dx -Rt)/dx)) /(grr*2.0);
        //float1 Trrr=(dot((Rr_dx -Rr)/dx,Rr) +dot(Rr,(Rr_dx -Rr)/dx)) /(grr*2.0);
//        float1 Trpp=(dot((Rp_dx -Rp)/dx,Rp) +dot(Rp,(Rp_dx -Rp)/dx)) /(grr*2.0);
        //float1 Tppp=(dot((Rp_dx -Rp)/dx,Rp) +dot(Rp,(Rp_dx -Rp)/dx)) /(gpp*2.0);
        
/*
Cuuu= <ruu;ru>  /<ru;ru>
Cvuu= <ruu;rv>  /<rv;rv>

Cuuv= <ruv;ru>  /<ru;ru>
Cvuv= <ruv;rv>  /<rv;rv>

Cuvv= <rvv;ru>  /<ru;ru>
Cvvv= <rvv;rv>  /<rv;rv>


ruu=Cuuu*ru + Cvuu*rv + L*normal
ruv=Cuuv*ru + Cvuv*rv + M*normal
rvv=Cuvv*ru + Cvvv*rv + N*normal
duu=Cuuu*ru + Cuuv*ru + Cuvv*ru    first row

        float1 Tttt=(dot((Rt_dx -Rt)/dx,Rt) +dot(Rt,(Rt_dx -Rt)/dx)) /(gtt*2.0);
        float1 Trtt=(dot((Rt_dx -Rt)/dx,Rt) +dot(Rt,(Rt_dx -Rt)/dx)) /(grr*2.0);
        float1 Trrr=(dot((Rr_dx -Rr)/dx,Rr) +dot(Rr,(Rr_dx -Rr)/dx)) /(grr*2.0);
        float1 Trpp=(dot((Rp_dx -Rp)/dx,Rp) +dot(Rp,(Rp_dx -Rp)/dx)) /(grr*2.0);
        float1 Tppp=(dot((Rp_dx -Rp)/dx,Rp) +dot(Rp,(Rp_dx -Rp)/dx)) /(gpp*2.0);


        float1 Tttt=(dot((wt_x-wt)/x,wt) +dot(wt,(wt_x-wt)/x)) /(gtt*2.0);
        float1 Trtt=(dot((wt_x-wt)/x,wt) +dot(wt,(wt_x-wt)/x)) /(grr*2.0);
        float1 Trrr=(dot((wr_x-wr)/x,wr) +dot(wr,(wr_x-wr)/x)) /(grr*2.0);
        float1 Trpp=(dot((wp_x-wp)/x,wp) +dot(wp,(wp_x-wp)/x)) /(grr*2.0);
        float1 Tppp=(dot((wp_x-wp)/x,wp) +dot(wp,(wp_x-wp)/x)) /(gpp*2.0);*/
        
        ddT  = -Ttrt*dT*dr*2;
        ddr  = -Trtt*dT*dT +Trrr*dr*dr  +Trpp*dp*dp;
        ddp  = -Tprp*dr*dp*2;


        r+=dr*dt;     
        p+=dp*dt;     
        t+=dt;         

        dr+=ddr*dt;
        dp+=ddp*dt;
        dT+=ddT*dt;

        T+=dt/dT;    
        

        pos.x = r*cos(p)*scale;
        pos.y = r*sin(p)*scale;

        
        pixel5(pos,0xff0000);
    }
}
void schwarzschild3()
{
    float1 r,p,t,T, dr,dp, ddr,ddT,ddp,  dT;
    
    c2  = c ;
    rs2 = rs;


    r  = sqrt(pos.x*pos.x + pos.y*pos.y);
    p  = atan2(pos.y,pos.x);
    t = 0.0;

    dr = (pos.x * vel.x + pos.y * vel.y) / r;
    dp = (pos.x * vel.y - pos.y * vel.x) / (r*r);
    dT=1.0;

    ddr=0;
    ddp=0;
    ddT=0;
    
    T = 0.0;
    


    for(int j = 0;j<n;j++)
    {
        if (r < rs2) return;

//ds^2=((r-2u)/r) dt^2 - (r/(r-2u)) dr^2 - (r^2) do^2 -(r^2 sin2 o) df^2

        vec3 Rt(sqrt((r-rs)/r)*c,0.0L,0.0L);//first derivative
        vec3 Rr(0.0,sqrt(r/(r-rs)),0.0);
        vec3 Rp(0.0,0.0,sqrt(r*r));

		float1 gtt=dot(Rt,Rt);//metric tensor
		float1 grr=dot(Rr,Rr);
		float1 gpp=dot(Rp,Rp);
		


        float1 dx=1e-2;//difference for second derivative
        float1 r2=r+dx;
        vec3 Rt_dx(sqrt((r2-rs)/r2)*c,0.0L,0.0L);//neightb
        vec3 Rr_dx(0.0,sqrt(r2/(r2-rs)),0.0);
        vec3 Rp_dx(0.0,0.0,sqrt(r2*r2));
        
        
		vec3 Rtr=(Rt_dx -Rt)/dx;//second derivative of t respect to r
		vec3 Rrr=(Rr_dx -Rr)/dx;
		vec3 Rpr=(Rp_dx -Rp)/dx;


        float1 Ttrt=dot(Rtr,Rt)/gtt;
        float1 Trtt=dot(Rtr,Rt)/grr;//dot(Rtt,Rr)/grr  direction of derivation and direction of component exchange
        float1 Trrr=dot(Rrr,Rr)/grr;//dot(Rrr,Rr)/grr
        float1 Trpp=dot(Rpr,Rp)/grr;//dot(Rpp,Rr)/grr
        float1 Tprp=dot(Rpr,Rp)/gpp;

        
        ddT  = -Ttrt*dT*dr*2;
        ddr  = -Trtt*dT*dT +Trrr*dr*dr  +Trpp*dp*dp;
        ddp  = -Tprp*dr*dp*2;


        r+=dr*dt;     
        p+=dp*dt;     
        t+=dt;         

        dr+=ddr*dt;
        dp+=ddp*dt;
        dT+=ddT*dt;

        T+=dt/dT;    
        

        pos.x = r*cos(p);
        pos.y = r*sin(p);

        pixel5(pos,0xffff00);
    }
}

void newton()
{
    for(int j = 0;j<n;j++)
    {
        pos += vel*dt;

        vec3 av = pos;
        float1 r=length(pos);
        av/=r;
        
        float1 a = -m*g/(r*r);
        vel += av*a*dt;

        pixel5(pos,0x0000ff);
    }
}
void refractold()
{
    float1 r,c1,c2, cx,cy,w;
    vec3 pos4,vel4;
    

      pos.z=0;
      vel.z=0;

    r=length(pos)-rs;
    w=2*m*g/(r*c*c*4);    
    c1=c*(1-w)/((1+w)); 

 
    pos4=pos;
    vel4=vel;
    vel4.z=sqrtl(c1*c1 - vel.x*vel.x - vel.y*vel.y);
    
    for(int j = 0;j<n;j++)
    {
        pos4 = pos4+vel4*dt;
        pos=pos4;
        pos.z=0;

        vec3 N = normalize(pos);
        vec3 T = normalize(vel4);          
        T=normalize(T-N*dot(N,T));

        r=length(pos)-rs;
        w=2*m*g/(r*c*c*4);
        c2=c*(1-w)/((1+w));

      
        cx=dot(T,vel4);
        cy=dot(N,vel4);
        float1 cy2=cy;

        cx=cx*c2*c2/(c1*c1);
        cy=(c2*c2-cx*cx);
        if(cy<0.0) cy=-sqrtl(-cy);
        else       cy= sqrtl(cy);
        if(cy2<0.0) cy=-cy;

        vel4=normalize(T*cx + N*cy)*c2;
        c1=c2;

        pixel5(pos,0x00ffff);
    }
}
void refract()
{
    float1 r,c1,c2, cx,cy,w;
    vec3 pos4,vel4;
    
rs*=0.52;//0.52

      pos.z=0;
      vel.z=0;

	float1 w2=1.0;
    r=length(pos)-rs*w2;
//    w=2*m*g/(r*c*c*4);        c1=c*(1-w)/((1+w)); 
	c1=c*(1.0-rs/r);//*0.5

 
    pos4=pos;
    vel4=vel;
    vel4.z=sqrtl(c1*c1 - vel.x*vel.x - vel.y*vel.y);//c1
    
    for(int j = 0;j<n;j++)
    {
        pos4 = pos4+vel4*dt;
        pos=pos4;
        pos.z=0;

        vec3 N = normalize(pos);
        vec3 T = normalize(vel4);          
        T=normalize(T-N*dot(N,T));

        r=length(pos)-rs*w2;
//        w=2*m*g/(r*c*c*4);        c2=c*(1-w)/((1+w));
		c2=c*(1.0-rs/r);//*0.5

      
        cx=dot(T,vel4);
        cy=dot(N,vel4);
        float1 cy2=cy;

        cx=cx*c2*c2/(c1*c1);
        cy=(c2*c2-cx*cx);
        if(cy<0.0) cy=-sqrtl(-cy);
        else       cy= sqrtl(cy);
        if(cy2<0.0) cy=-cy;

        vel4=normalize(T*cx + N*cy)*c2;
        c1=c2;

        pixel5(pos,0x00ffff);
    }
}
vec3 getsurfacepoint(float1 u,float1 v)
{
	vec3 v1;
	float1 r=sqrtl(u*u+v*v);
	v1.x=u;
	v1.y=v;

	v1.z=2*sqrtl(rs*(r*2-rs));	//Flamm paraboloid
	v1.z*=c/length(vel);
	
	return v1;
}
float1 h=0;

vec3 getsurfacenormal(vec3 p0)
{
	float1 u=p0.x;
	float1 v=p0.y;
	float1 d=1e-5;
	
	vec3 p1=getsurfacepoint(u,v  );
	vec3 p2=getsurfacepoint(u+d,v);
	vec3 p3=getsurfacepoint(u,v+d);
	
	return normalize(cross(p2-p1,p3-p1));
}
void flammsurface()
{
	for(int i=0;i<n;i++)
	{
		vec3 normal=getsurfacenormal(pos);
		pos+=vel*dt;
		vel.z-=h*dt*(1.0-sqrt(fabs(normal.z)));
		vel-=normal*dot(normal,vel);

        pixel5(pos,0x00ff00);//green
	}
}

void fluidgrav()
{
	float f0,mf,r2,dr,dt2,v11,v22,f2_doppler,f1_doppler,l1,l2,p1,p2,dp,dv,a,h=6.626e-34;
   	
   	f0=1e15; //some frequency 
	mf=f0*h/(c*c); //f=E*h  equivalent mass
//	dr=rs*30.0;//0.02
	dr=rs*0.001;//0.02  0.0001-0.1
	dt2=2.0*dr/c;
		
    for(int j = 0;j<n;j++)
    {
        pos += vel*dt;

        vec3 av = pos;
        float1 r=length(pos);
        av/=r;

// rs*3????????
#if 0
		r2=r+dr;	v11=c*(rs/(2.0*r2 -rs*3.0));//2 value of velocity field near to test mass
		r2=r-dr;	v22=c*(rs/(2.0*r2 -rs*3.0));
#else
//		r2=r+dr;	v11=c/sqrt(1.0 -rs/r2)-c;//original
//		r2=r-dr;	v22=c/sqrt(1.0 -rs/r2)-c;
		r2=r+dr;	v11=c/sqrt(1.0 -rs/(r2-rs))-c;//rs origo, GOOD!
		r2=r-dr;	v22=c/sqrt(1.0 -rs/(r2-rs))-c;
#endif

#if 0
		f1_doppler=f0*(c+v11)/c;//normal doppler
		f2_doppler=f0*(c+v22)/c;

		l1=c/f1_doppler;  //wavelength
		l2=c/f2_doppler; 
		p1=h/l1;//moment         p1=h/(c/f1_doppler)=f1_doppler*h/c;
		p2=h/l2;
	
		dp=(p2-p1); // delta moment
		dv=dp/mf; //delta velocity
		dv=(p2-p1)/mf;
		a=-dv/dt2;//acceleration
#else
		dv=v22-v11;
		a=-dv/dt2;//acceleration
#endif
        vel += av*a*dt;//        float1 a = -m*g/(r*r);

        pixel5(pos,0xff00ff);
    }
}

void reset()
{

     g = 6.674e-11;
     c = 2.997e8;
     c2 =1.0;
     m = 1.98e30;
float1 s=1.1;

    rs = 2.0*m*g/(c*c);
    r1 = rs*1500.0;    dt = 2e-5*s;
//    r1 = rs*150.0;    dt = 2e-7*s;
//    r1 = rs*60.0*s;    dt = 2e-8*s;
//    r1 = rs*22.0*s;    dt =5e-9*s;
//    r1 = rs*15.0*s;    dt =2e-9*s;
    v1 = sqrt(m*g/r1)*1.2;
    
    n = 2000000*3;
    scale2 = 300.0/r1;
   
    pos = vec3(r1,(float1)0.0,(float1)0.0);
    vel = vec3(v1*(float1)0.3,v1*(float1)0.7,(float1)0.0);
}
int main()
{
	init_system();   

	qclear();
	
    reset();        newton();   
    reset();        schwarzschild();
    reset();        fluidgrav();
    
//    reset();        schwarzschild2();
//    reset();        schwarzschild3();
//    reset();		refract();
//    reset();		flammsurface();
    
    qflush();
    getchar();

    return 0;

}


