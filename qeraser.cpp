
//g++ x.cpp -lX11 -O3 

#define __LINUX__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include <pthread.h>
#include <X11/Xlib.h>
#include <X11/Xatom.h>                                                             
#include <X11/keysym.h>                                                                   
#include <fcntl.h>



Display *dpy;
Window win;
GC gc;


#define iradian (M_PI/180.0)
#define frnd(n) (frand()*(float)n)



float sqr(float n)
{
	return n*n;
}
float frand() //rand  0.0-1.0
{
    return (float)(rand()%10000)/10000.0;
}





//complex or ordinary 2d vector
#if 0

struct vec2
{
    float x,y;
    
    vec2() {x=0;y=0;}
    vec2(float a,float b) {x=a;y=b;}
    vec2 operator -(vec2 a)
    {
    	vec2 b;
    	b.x=x - a.x;
    	b.y=y - a.y;
    	return b;
    }
    vec2 operator *(float c)
    {
    	vec2 u;
    	u.x=x*c;
    	u.y=y*c;
    	return u;
    }
};
float dot(vec2 v1,vec2 v2)
{
    return v1.x*v2.x + v1.y*v2.y;
}
void add_amp(vec2 *v,float phase,float ax)
{
    v->x += sin(phase)*ax;
    v->y += cos(phase)*ax;
};
void add_amp(vec2 *v,float phase,float ax,float ay)
{
    v->x += sin(phase)*ax;
    v->y += cos(phase)*ay;
};
void add_polarizer2(vec2 *v,float phase)
{
	vec2 pol;

	add_amp(&pol,phase,1.0);
	float amp=dot(*v,pol);      //v'=pol><pol|v>
    v->x = pol.x*amp;    
    v->y = pol.y*amp;
};
float probability(vec2 v1,vec2 v2)
{
    return dot(v1,v2); // length(v1)^2
}
float length(vec2 v1)
{
    return sqrt(dot(v1,v1));
}
void rot(vec2 &v1)
{
	float a=v1.x;
	v1.x=v1.y;
	v1.y=-a;
}
void normalize(vec2 v1)
{
	float e=sqrt(dot(v1,v1));
	
	v1.x/=e;
	v1.y/=e;
}


//linear to circular
void add_quarterwaveplatel2c(vec2 *v,float phase,float dist,float qwp_angle,float ax)
{
    vec2 axis_fast,axis_slow,input_wave;
    
    add_amp(&axis_fast,qwp_angle,1.0);
    add_amp(&axis_slow,qwp_angle+90.0*iradian,1.0);
    add_amp(&input_wave,phase,1.0);

    float amp_fast=dot(axis_fast,input_wave);
    float amp_slow=dot(axis_slow,input_wave);
    
    amp_fast=cos(dist)*amp_fast;
    amp_slow=sin(dist)*amp_slow;//-+90 phase shift
    
    v->x += (axis_fast.x*amp_fast + axis_slow.x*amp_slow)*ax;
    v->y += (axis_fast.y*amp_fast + axis_slow.y*amp_slow)*ax;
}


//__________________________________________________________________
#else
#define vec2 complex

struct complex
{
	union
	{
    	float real;
    	float x;
    };
	union
	{
    	float img;
    	float y;
    };

    complex() {real=img=0;};
    complex(float r,float i) {real=r;img=i;};
    void set(float r,float i) {real=r;img=i;};
    void set(complex c) {real=c.real;img=c.img;};

    void conjugate() {img*=-1.0;}
    complex operator +(complex c)
    {
        complex e;
        e.real = this->real + c.real;
        e.img = this->img + c.img;
        return e;
    }
    complex operator -(complex c)
    {
        complex e;
        e.real = this->real - c.real;
        e.img = this->img - c.img;
        return e;
    }
    complex operator *(complex c)
    {
        complex e;
        e.real = this->real*c.real - this->img*c.img;
        e.img = this->img*c.real + this->real*c.img;
        return e;
    }
    complex operator /(complex c)
    {
        complex e;
        float denominator=c.real*c.real + c.img*c.img;
        e.real = (this->real * c.real + this->img * c.img)/denominator;
        e.img = (this->img * c.real - this->real * c.img)/denominator;
        return e;
    }
    complex operator *(float c)
    {
        complex e;
        e.real = this->real*c;
        e.img = this->img*c;
        return e;
    }
};
void add_amp(complex *v,float phase,float ax)
{
    *v = *v + complex(sin(phase)*ax,cos(phase)*ax);
};
void add_polarizer2(complex *v,float phase)
{
	complex pol,pol2;
	add_amp(&pol,phase,1.0);
	pol2=pol;
	pol2.conjugate();
	complex amp=(*v)*pol2; amp.img=0;      //v'=pol><pol|v>
    *v = pol*amp;    
};
float probability(complex c1,complex c2)
{
    c1.conjugate();// <c1|
    complex P=c1*c2; //<c1|c2>

    return P.real;
}
float length(complex c1)
{
	complex c2=c1;
    c1.conjugate();// <c1|
    complex P=c1*c2; //<c1|c2>

    return sqrt(P.real);
}
void add_quarterwaveplatel2c(vec2 *v,float phase,float dist,float qwp_angle,float ax)
{
    complex axis_fast,axis_slow,input_wave;
    
    add_amp(&axis_fast,qwp_angle,1.0);
    add_amp(&axis_slow,qwp_angle+90.0*iradian,1.0);
    add_amp(&input_wave,phase,1.0);

	input_wave.conjugate();
    complex amp_fast=axis_fast*input_wave; amp_fast.img=0;
    complex amp_slow=axis_slow*input_wave; amp_slow.img=0;
    
    amp_fast=amp_fast*complex(cos(dist)*ax,0.0);
    amp_slow=amp_slow*complex(sin(dist)*ax,0.0);//-+90 phase shift
    
    *v=*v+axis_fast*amp_fast + axis_slow*amp_slow;
};

#endif






int wx=1100;
int wy=650;

XEvent event;               
int mbut=0,mx,my,dmx,dmy,key=0;




void pixel(int x,int y,int col)
{
    XSetForeground(dpy,gc,col);
    XDrawPoint(dpy, win, gc, x,y);
}
void pixel(int x,int y,int r,int g,int b)
{
    XSetForeground(dpy,gc,(r<<16)+(g<<8)+b);
    XDrawPoint(dpy, win, gc, x,y);
}
void line(int x1,int y1,int x2,int y2,int col)
{
    XSetForeground(dpy,gc,col);
    XDrawLine(dpy, win, gc, x1,y1,x2,y2);
}
void line(int x1,int y1,int x2,int y2,int r,int g,int b)
{
    XSetForeground(dpy,gc,(r<<16)+(g<<8)+b);
    XDrawLine(dpy, win, gc, x1,y1,x2,y2);
}



int init_system(int _wx=0,int _wy=0)
{

	if(_wx!=0) wx=_wx;
	if(_wy!=0) wy=_wy;

	dpy = XOpenDisplay(0);
	if(dpy==0) printf("Do not use in root ! \n");
	win = XCreateSimpleWindow(dpy, DefaultRootWindow(dpy), 0,0, wx, wy, 0,0,0);

	XSelectInput(dpy, win, StructureNotifyMask| ButtonPressMask| ButtonReleaseMask| KeyPressMask |PointerMotionMask);
	XMapWindow(dpy,win);
	gc = XCreateGC(dpy, win, 0, 0);

#if 0
    const char *pattern="*15*";//font manager   12,14,18
    int maxnames=10;
    int actual_count_return=0;
	char **fontnames=XListFonts(dpy, pattern, maxnames, &actual_count_return);
	for(int i=0;i<actual_count_return;i++)	printf("%s %d \n",fontnames[i],actual_count_return);
		
	XFontStruct *fontinfo=0; 
	fontinfo = XLoadQueryFont(dpy,fontnames[0]); 
	XSetFont(dpy,gc,fontinfo->fid); 
#endif



	for(;;) { XEvent e; XNextEvent(dpy, &e); if (e.type == MapNotify) break; }
	XMapWindow(dpy, win);


	return 0;
}

int keyboard()
{
//	mbut=0;
	int u;
	
	dmx=0;
	dmy=0;

        while (XPending(dpy) > 0)        
        {                                    
            XNextEvent(dpy, &event);     
//            XPeekEvent(disp, &event);     
            switch (event.type)              
            {                                
                case Expose:                 
                    break;                       
                case ButtonPress:
                	mx=event.xbutton.x;
	                my=event.xbutton.y;
//	                printf("%d %d \n",xm,ym);
					mbut=1;
                    break;
                case ButtonRelease:
                	mx=event.xbutton.x;
	                my=event.xbutton.y;
					mbut=0;
                    break;
                case MotionNotify:
                	dmx=mx-event.xmotion.x;               
                	dmy=my-(event.xmotion.y);
                	mx=event.xmotion.x;               
                	my=event.xmotion.y;
                	
                	if(dmx ||dmy) return 1;
                    break;
                case KeyPress:
                	u=XLookupKeysym(&event.xkey, 0);
//	                printf("%d \n",u);

					key=u;
				
                    if (u == XK_Escape)
                    {
                    	exit(3);
                    }
                    break;
                case ClientMessage:
//                    if (strcmp(XGetAtomName(disp, event.xclient.message_type),
  //                             "WM_PROTOCOLS") == 0)
                    {
                    }
                    break;
                default:
                    break;
            }
        }
        return 0;
}



float cos2(float a) 
{
	a=cos(a);
	return a*a;
}

float sin2(float a) 
{
	a=sin(a);
	return a*a;
}



int dcqe()
{
	double scale=250.0;


    for(int ds_x=0;ds_x<200;ds_x++)//Ds position  -+4mm    
    {
        double P01=0,P02=0,P03=0,P04=0,
	        Pd0=0,Pd1=0,Pd2=0,Pd3=0,Pd4=0;
        
        double step=M_PI*2.0/300.0;
        for(double p=0;p<M_PI*2.0;p+=step)
        {
            vec2 amp_d0,amp_d1,amp_d2,amp_d3,amp_d4;
#if 0
            double photon_red=p;
            double photon_blue=photon_red+M_PI/2;
#else            
            double photon_red=M_PI*2*frand();
            double photon_blue=M_PI*2*frand();
#endif
            

            double ds_distance=850.0;
            double dp_distance=980.0;
            double wavelength=702.2e-6;
            double k1=2.0*M_PI/wavelength;//wave number
            double wavelength2=452.2e-6;
            double k2=2.0*M_PI/wavelength2;//wave number
            double hole_dist=0.3;
            double hole_wide=0.2;//0.2mm  = 200 micrometer wide   
            double ds_pos=4.0*(double)(ds_x-100)/100.0;//+-4mm     D0 position
            int slit_wide=10;
			
            for(int w=0;w<slit_wide;w++)
            {
                double hole1x=hole_dist/2.0 + hole_wide*(double)w/slit_wide;//hole A
                double hole2x=-hole_dist/2.0 - hole_wide*(double)w/slit_wide;//hole B
                double dist1=sqrt(sqr(ds_pos - hole1x) + sqr(ds_distance));
                double dist2=sqrt(sqr(ds_pos - hole2x) + sqr(ds_distance));    

                add_amp(&amp_d0,photon_red +dist1*k1 ,0.5);//2 slit  DS/D0 side
                add_amp(&amp_d0,photon_blue +dist2*k2 ,0.5);
            }
            double normalization_factor=1.0/(double)slit_wide;
            amp_d0=amp_d0*normalization_factor;

			// DP/D1D2D3D4 side
           	add_amp(&amp_d1,photon_red +M_PI+M_PI/2,0.5);//1 reflection + 1 transmit
           	add_amp(&amp_d1,photon_blue+M_PI*2 ,0.5);//2  reflection

            add_amp(&amp_d2,photon_blue+M_PI+M_PI/2 ,0.5);//1  reflection+ 1 transmit
            add_amp(&amp_d2,photon_red +M_PI*2,0.5);//2 reflection
            
            add_amp(&amp_d3,photon_blue+M_PI ,1.0);//1 reflection
            
            add_amp(&amp_d4,photon_red +M_PI,1.0);//1 reflection

            P01+=probability(amp_d1,amp_d1)*probability(amp_d0,amp_d0)*step;//common detections
            P02+=probability(amp_d2,amp_d2)*probability(amp_d0,amp_d0)*step;
            P03+=probability(amp_d3,amp_d3)*probability(amp_d0,amp_d0)*step;
            P04+=probability(amp_d4,amp_d4)*probability(amp_d0,amp_d0)*step;

            Pd0+=probability(amp_d0,amp_d0)*step;// detectors alone
            Pd1+=probability(amp_d1,amp_d1)*step;
            Pd2+=probability(amp_d2,amp_d2)*step;
            Pd3+=probability(amp_d3,amp_d3)*step;
            Pd4+=probability(amp_d4,amp_d4)*step;
        }
        double normalization_factor2=(M_PI*2.0);
        P01/=normalization_factor2;      pixel(ds_x*2,200-P01*scale,0xffff00);//d0-d1234 common detection yellow
        P02/=normalization_factor2;      pixel(ds_x*2,320-P02*scale,0xffff00);
        P03/=normalization_factor2;      pixel(ds_x*2,440-P03*scale,0xffff00);
        P04/=normalization_factor2;      pixel(ds_x*2,550-P04*scale,0xffff00);
        
        Pd0/=normalization_factor2;     pixel(ds_x*2+400,550-Pd0*scale,0xff0000);//d0 alone red
        
        Pd1/=normalization_factor2;     pixel(ds_x*2+400,200-Pd1*scale,0x00ff);//d1234 alone  blue
        Pd2/=normalization_factor2;     pixel(ds_x*2+400,320-Pd2*scale,0x00ff);
        Pd3/=normalization_factor2;     pixel(ds_x*2+400,440-Pd3*scale,0x00ff);
        Pd4/=normalization_factor2;     pixel(ds_x*2+400,550-Pd4*scale,0x00ff);
    }

    XFlush(dpy);
    getchar();
        
    return 0;
}


int main()
{
	init_system();
	dcqe();
}

