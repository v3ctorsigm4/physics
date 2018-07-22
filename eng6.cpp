//g++ x.cpp -lX11 -O3 

#define __LINUX__

#ifndef __USE_X11
#define __USE_GLX
#endif
//#define  _USEMSAA



#ifdef __USE_GLX
#define GL_GLEXT_PROTOTYPES

#include <GL/glx.h>
#include <GL/gl.h>
#include <unistd.h>
#include <iostream>
 
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

 
#define GLX_CONTEXT_MAJOR_VERSION_ARB       0x2091
#define GLX_CONTEXT_MINOR_VERSION_ARB       0x2092
typedef GLXContext (*glXCreateContextAttribsARBProc)(Display*, GLXFBConfig, GLXContext, Bool, const int*);

GLXContext ctx=0;
#endif



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include <pthread.h>
#include <X11/Xlib.h>
#include <X11/Xatom.h>                                                             
#include <X11/keysym.h>                                                                   
#include <fcntl.h>


#include <string>
#include <vector>
#include <map>

using namespace std;

Display *display5;
Window win;
GC gc;


float delta_e=0.1;
int needarrowhead=1;


//#define float1 float
#define float1 double


#define iradian (M_PI/180.0)
#define frnd(n) (frand()*(float)n)


#include "complex.cpp"


float1 fsgn(float1 n)
{
	if(n<0.0) return -1.0;
	if(n>0.0) return 1.0;
	return 0.0;
}

//#define DIVLIMIT 1e-21
#define DIVLIMIT 1e-11

inline void chklim(float1 &n)
{
	if(n==0.0) n=DIVLIMIT;
	else	if(fabs(n)<DIVLIMIT) n=DIVLIMIT*fsgn(n);
}

struct vec4
{
	float1 x,y,z,w;

vec4 operator - (vec4 v) {vec4 v_;v_.x=x-v.x;v_.y=y-v.y;v_.z=z-v.z;v_.w=w-v.w;return v_;};
vec4 operator / (float1 s) {chklim(s);vec4 v_;v_.x=x/s;v_.y=y/s;v_.z=z/s;v_.w=w/s;return v_;};
vec4 operator - () {vec4 v_;v_.x=-x;v_.y=-y;v_.z=-z;v_.w=-w;return v_;};
};

struct vec3
{
float1 x,y,z;

vec3() {x=0;y=0;z=0;};
vec3(float1 x_,float1 y_,float1 z_) {x=x_;y=y_;z=z_;};

vec3 operator + (vec3 v) {vec3 v_;v_.x=x+v.x;v_.y=y+v.y;v_.z=z+v.z;return v_;};
vec3 operator - (vec3 v) {vec3 v_;v_.x=x-v.x;v_.y=y-v.y;v_.z=z-v.z;return v_;};
vec3 operator * (vec3 v) {vec3 v_;v_.x=x*v.x;v_.y=y*v.y;v_.z=z*v.z;return v_;};
vec3 operator / (vec3 v) {chklim(v.x);chklim(v.y);chklim(v.z);vec3 v_;v_.x=x/v.x;v_.y=y/v.y;v_.z=z/v.z;return v_;};

vec3 operator += (vec3 v) {x+=v.x;y+=v.y;z+=v.z;return *this;};
vec3 operator -= (vec3 v) {x-=v.x;y-=v.y;z-=v.z;return *this;};
vec3 operator *= (vec3 v) {x*=v.x;y*=v.y;z*=v.z;return *this;};
vec3 operator /= (vec3 v) {chklim(v.x);chklim(v.y);chklim(v.z);x/=v.x;y/=v.y;z/=v.z;return *this;};

vec3 operator + (float1 s) {vec3 v_;v_.x=x+s;v_.y=y+s;v_.z=z+s;return v_;};
vec3 operator - (float1 s) {vec3 v_;v_.x=x-s;v_.y=y-s;v_.z=z-s;return v_;};
vec3 operator * (float1 s) {vec3 v_;v_.x=x*s;v_.y=y*s;v_.z=z*s;return v_;};
vec3 operator / (float1 s) {chklim(s);vec3 v_;v_.x=x/s;v_.y=y/s;v_.z=z/s;return v_;};

vec3 operator += (float1 s) {x+=s;y+=s;z+=s;return *this;};
vec3 operator -= (float1 s) {x-=s;y-=s;z-=s;return *this;};
vec3 operator *= (float1 s) {x*=s;y*=s;z*=s;return *this;};
vec3 operator /= (float1 s) {chklim(s);x/=s;y/=s;z/=s;return *this;};

vec3 operator - () {vec3 v_;v_.x=-x;v_.y=-y;v_.z=-z;return v_;};
};

float1 dot(vec3 v1,vec3 v2) { return (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);}
float1 length(vec3 v1) 
{
	float1 t1=dot(v1,v1);
	if(t1<DIVLIMIT) return t1; 
	return sqrt(t1);
}
float1 length4(vec3 v1,float1 w) 
{
	float1 t1=dot(v1,v1) + w*w;
	if(t1<DIVLIMIT) return t1; 
	return sqrt(t1);
}
vec3 normalize(vec3 v1) 
{ 
	float1 t1=length(v1);
	chklim(t1);
	return v1/t1;
}
vec3 cross(vec3 va,vec3 vb)
{
vec3 v_;

v_.x=(va.y*vb.z)-(vb.y*va.z);
v_.y=(va.z*vb.x)-(vb.z*va.x);
v_.z=(va.x*vb.y)-(vb.x*va.y);

return v_;
};

void cvariant(vec3 &v1) 	{	v1=v1/dot(v1,v1);}; 

inline float1 clamp(float1 n,float1 l1,float1 l2)
{
	if(n<l1) n=l1;
	if(n>l2) n=l2;
	
	return n;
}
inline vec3 clamp(vec3 n,float1 l1,float1 l2)
{
	n.x=clamp(n.x,l1,l2);	
	n.y=clamp(n.y,l1,l2);	
	n.z=clamp(n.z,l1,l2);	
	
	return n;
}
inline float1 saturate(float1 n)
{
	if(n<0.0) n=0.0;
	if(n>1.0) n=1.0;
	
	return n;
}
inline vec3 saturate(vec3 v1)
{
	v1.x=saturate(v1.x);
	v1.y=saturate(v1.y);
	v1.z=saturate(v1.z);
	
	return v1;
}
inline float1 lerp(float1 v1,float1 v2,float1 t)
{
	return v1+(v2-v1)*t;
}
inline vec3 lerp(vec3 v1,vec3 v2,float1 t)
{
	v1.x=lerp(v1.x,v2.x,t);
	v1.y=lerp(v1.y,v2.y,t);
	v1.z=lerp(v1.z,v2.z,t);
	
	return v1;
}
inline vec3 lerp(vec3 v1,vec3 v2,vec3 t)
{
	v1.x=lerp(v1.x,v2.x,t.x);
	v1.y=lerp(v1.y,v2.y,t.y);
	v1.z=lerp(v1.z,v2.z,t.z);
	
	return v1;
}

float1 sqr(float1 n)
{
	return n*n;
}
float1 frand()
{
    return (float1)(rand()%10000)/10000.0;
}




struct axis
{
	vec3 x;
	vec3 y;
	vec3 z;
	
	axis()
	{
		x=vec3(1,0,0);
		y=vec3(0,1,0);
		z=vec3(0,0,1);
	}
};


int wx=1100;
int wy=650;

XEvent event;               
int mbut=0,mx,my,dmx,dmy,key=0;
int viewmode=1;
vec3 camx,camy,camz,eye,look,up;

vec3 coltab[6];
float1 asp=1.0;





#ifdef __USE_GLX

int lineCC=0;


void qclear()
{
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glLoadIdentity();
}

void qflush()
{
	if(lineCC>0)
	{
		glEnd();
		lineCC=0;
	}
    glXSwapBuffers (display5, win);
}


void pixel(int x,int y,int r,int g,int b)
{
	if(lineCC==0) glBegin(GL_LINES);

	float1 r1=(float1)r/255.0;
	float1 g1=(float1)g/255.0;
	float1 b1=(float1)b/255.0;
	
	float1 x1b=asp*(float1)x/1000.0-wx/2;
	float1 y1b=    (float1)y/1000.0-wy/2;
	float1 x2b=x1b+0.001;
	float1 y2b=y1b;
		
	glColor3f(r1,g1,b1);	glVertex3f(x1b,y1b,0.0);
	glColor3f(r1,g1,b1);	glVertex3f(x2b,y2b,0.0);
	lineCC+=2;

	if(lineCC>16000)
	{
		glEnd();
		lineCC=0;
	}
}
void pixel(int x,int y,int col)
{
	int r=(col>>16);
	int g=(col>>8)&255;
	int b=(col)&255;
	pixel(x,y,r,g,b);
}
void line(int x1,int y1,int x2,int y2,int r,int g,int b)
{
	if(lineCC==0) glBegin(GL_LINES);

	float1 r1=(float1)r/255.0;
	float1 g1=(float1)g/255.0;
	float1 b1=(float1)b/255.0;
	
	float1 x1b=(float1)x1/1000.0;
	float1 x2b=(float1)x2/1000.0;
	float1 y1b=(float1)y1/1000.0;
	float1 y2b=(float1)y2/1000.0;
		
	glColor3f(r1,g1,b1);	glVertex3f(x1b,y1b,0.0);
	glColor3f(r1,g1,b1);	glVertex3f(x2b,y2b,0.0);
	lineCC+=2;

	if(lineCC>16000)
	{
		glEnd();
		lineCC=0;
	}
}
void line(int x1,int y1,int x2,int y2,int col)
{
	int r=(col>>16);
	int g=(col>>8)&255;
	int b=(col)&255;
	line(x1,y1,x2,y2,r,g,b);
}



void init_system(int _wx=0,int _wy=0)
{
	
	if(_wx!=0) wx=_wx;
	if(_wy!=0) wy=_wy;

	asp=(float1)wx/(float1)wy;

        display5 = XOpenDisplay(0);
        if(display5==0) printf("ERROR XOpenDisplay(0) \n");
  	
  
        glXCreateContextAttribsARBProc glXCreateContextAttribsARB = NULL;
     
        const char *extensions = glXQueryExtensionsString(display5, DefaultScreen(display5));
        std::cout << extensions << std::endl;
     
#ifdef  _USEMSAA
    static int visual_attribs[] =
    {
       GLX_X_RENDERABLE    , True,
      GLX_DRAWABLE_TYPE   , GLX_WINDOW_BIT,
      GLX_RENDER_TYPE     , GLX_RGBA_BIT,
      GLX_X_VISUAL_TYPE   , GLX_TRUE_COLOR,
      GLX_RED_SIZE        , 8,
      GLX_GREEN_SIZE      , 8,
      GLX_BLUE_SIZE       , 8,
      GLX_ALPHA_SIZE      , 8,
      GLX_DEPTH_SIZE      , 24,
      GLX_STENCIL_SIZE    , 8,
      GLX_DOUBLEBUFFER    , True,//True,
      GLX_SAMPLE_BUFFERS  , 1,            // <-- MSAA
      GLX_SAMPLES         , 2,            // <-- MSAA
      None
	};
#else	
    static int visual_attribs[] =
    {
       GLX_X_RENDERABLE    , True,
        GLX_RENDER_TYPE, GLX_RGBA_BIT,
        GLX_DRAWABLE_TYPE, GLX_WINDOW_BIT,
      GLX_X_VISUAL_TYPE   , GLX_TRUE_COLOR,
        GLX_DOUBLEBUFFER, true,
        GLX_RED_SIZE, 8,
        GLX_GREEN_SIZE, 8,
        GLX_BLUE_SIZE, 8,
      GLX_DEPTH_SIZE      , 24,
      GLX_STENCIL_SIZE    , 8,
        None
     };     
#endif


        std::cout << "Getting framebuffer config" << std::endl;
        int fbcount;
        GLXFBConfig *fbc = glXChooseFBConfig(display5, DefaultScreen(display5), visual_attribs, &fbcount);
        if (!fbc)
        {
            std::cout << "Failed to retrieve a framebuffer config" << std::endl;
            return ;
        }
     
        std::cout << "Getting XVisualInfo" << std::endl;
        XVisualInfo *vi = glXGetVisualFromFBConfig(display5, fbc[0]);
/*		for(int i=0;i<fbcount;i++)
		{
	        XVisualInfo *vi2 = glXGetVisualFromFBConfig(display5, fbc[i]);
				printf("%d  %d %p\n",fbcount,(int)vi2->screen,vi2->visual);
		}*/
     
        XSetWindowAttributes swa;
        std::cout << "Creating colormap" << std::endl;
        swa.colormap = XCreateColormap(display5, RootWindow(display5, vi->screen), vi->visual, AllocNone);
//  swa.background_pixmap = None ;//new
//  swa.background_pixel  = 0    ;
          swa.border_pixel = 0;
        swa.event_mask = StructureNotifyMask;
   
        std::cout << "Creating window" << std::endl;//yy-90
    win = XCreateWindow(display5, RootWindow(display5, vi->screen), 0, 0, wx, wy, 0, vi->depth, InputOutput, vi->visual, CWBorderPixel|CWColormap|CWEventMask, &swa);        
        if (!win)
        {
            std::cout << "Failed to create window." << std::endl;
            return ;
        }
     
        std::cout << "Mapping window" << std::endl;
        XMapWindow(display5, win);


        // Create an oldstyle context first, to get the correct function pointer for glXCreateContextAttribsARB
        GLXContext ctx_old = glXCreateContext(display5, vi, 0, GL_TRUE);//TRUE
        glXCreateContextAttribsARB =  (glXCreateContextAttribsARBProc)glXGetProcAddress((const GLubyte*)"glXCreateContextAttribsARB");
        glXMakeCurrent(display5, 0, 0);
        glXDestroyContext(display5, ctx_old);

        if (glXCreateContextAttribsARB == NULL)
        {
            std::cout << "glXCreateContextAttribsARB entry point not found. Aborting." << std::endl;
            return ;
        }
     
        static int context_attribs[] =
        {
            GLX_CONTEXT_MAJOR_VERSION_ARB, 2,//3
            GLX_CONTEXT_MINOR_VERSION_ARB, 0,//0
            None
        };

        std::cout << "Creating context" << std::endl;
        ctx = glXCreateContextAttribsARB(display5, fbc[0], NULL, true, context_attribs);
        if (!ctx)
        {
            std::cout << "Failed to create GL3 context." << std::endl;
            return ;
        }

        std::cout << "Making context current" << std::endl;
        glXMakeCurrent(display5, win, ctx);

	XSelectInput(display5, win, StructureNotifyMask| 
			ButtonPressMask| ButtonReleaseMask| 
			KeyPressMask |KeyReleaseMask |
			PointerMotionMask);


		

        std::cout << "Making context OK" << std::endl;

//            glClearColor (0.0, 1.0, 1.0, 1.0);            glClear (GL_COLOR_BUFFER_BIT);           // glXSwapBuffers (display5, win);
       
//             sleep(1);           

 
#ifdef  _USEMSAA
 glEnable( GL_MULTISAMPLE ); 
 glEnable(GL_POLYGON_SMOOTH); 
 glEnable(GL_LINE_SMOOTH);
 glEnable(GL_POINT_SMOOTH) ;
#endif

   	glViewport(0,0, wx, wy);

	//ERR(glDisable(GL_DEPTH));	
	glDisable (GL_DEPTH_TEST);
	glDisable (GL_BLEND);
	glDisable(GL_TEXTURE_2D);
	

        std::cout << "init OK" << std::endl;

}


#else

void qclear()
{
	XClearWindow(display5,win);
}

void qflush()
{
	XFlush(display5);
}

void pixel(int x,int y,int col)
{
    XSetForeground(display5,gc,col);
    XDrawPoint(display5, win, gc, x,y);
}
void pixel(int x,int y,int r,int g,int b)
{
    XSetForeground(display5,gc,(r<<16)+(g<<8)+b);
    XDrawPoint(display5, win, gc, x,y);
}
void line(int x1,int y1,int x2,int y2,int col)
{
    XSetForeground(display5,gc,col);
    XDrawLine(display5, win, gc, x1,y1,x2,y2);
}
void line(int x1,int y1,int x2,int y2,int r,int g,int b)
{
    XSetForeground(display5,gc,(r<<16)+(g<<8)+b);
    XDrawLine(display5, win, gc, x1,y1,x2,y2);
}


void init_system(int _wx=0,int _wy=0)
{
	
	if(_wx!=0) wx=_wx;
	if(_wy!=0) wy=_wy;

	asp=(float1)wx/(float1)wy;
	
	display5 = XOpenDisplay(0);
	if(display5==0) printf("Do not use in root ! \n");
	win = XCreateSimpleWindow(display5, DefaultRootWindow(display5), 0,0, wx, wy, 0,0,0);

	XSelectInput(display5, win, StructureNotifyMask| ButtonPressMask| ButtonReleaseMask| KeyPressMask |PointerMotionMask);
	XMapWindow(display5,win);
	gc = XCreateGC(display5, win, 0, 0);

#if 0
    const char *pattern="*15*";//font manager   12,14,18
    int maxnames=10;
    int actual_count_return=0;
	char **fontnames=XListFonts(display5, pattern, maxnames, &actual_count_return);
	for(int i=0;i<actual_count_return;i++)	printf("%s %d \n",fontnames[i],actual_count_return);
		
	XFontStruct *fontinfo=0; 
	fontinfo = XLoadQueryFont(display5,fontnames[0]); 
	XSetFont(display5,gc,fontinfo->fid); 
#endif


	for(;;) { XEvent e; XNextEvent(display5, &e); if (e.type == MapNotify) break; }
	XMapWindow(display5, win);


	return ;
}
#endif


int keyboard()
{
//	mbut=0;
	int u;
	
	dmx=0;
	dmy=0;

        while (XPending(display5) > 0)        
        {                                    
            XNextEvent(display5, &event);     
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
					if(u==49) viewmode=1;
					if(u==50) viewmode=2;
					if(u==51) viewmode=3;
					if(u==52) viewmode=4;
					if(u==53) viewmode=5;
					if(u==54) viewmode=6;
					if(u==55) viewmode=7;
					if(u==56) viewmode=8;
				
                    if (u == XK_Escape)
                    {
	                     XDestroyWindow(display5,win);                    	
#ifdef __USE_GLX
					    ctx = glXGetCurrentContext();
					    glXMakeCurrent(display5, 0, 0);
    					glXDestroyContext(display5, ctx);
#endif
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








int transf3d2(vec3 v1,int *xc,int *yc)
{
	vec3 w1=v1-eye;
	float1 zoom=750;//650
	
	float1 x=dot(w1,camx);
	float1 y=dot(w1,camy);
	float1 z=dot(w1,camz);
#ifdef __USE_GLX
	xc[0]= x*zoom/z;//perspective projection
	yc[0]= y*zoom*asp/z;
#else
	xc[0]=wx/2+ x*zoom/z;//perspective projection
	yc[0]=wy/2- y*zoom/z;
#endif
	if(z<0.0) return 0;
	
	return 1;
}
void line3d(vec3 p1,vec3 p2,int color)
{
	int x1,y1,x2,y2;

	if(transf3d2(p1,&x1,&y1))
	if(transf3d2(p2,&x2,&y2))
		line(x1,y1,x2,y2,color);
}
void line2d(vec3 p1,vec3 p2,int color)
{
	int x3=wx/2,y3=wy/2;

	p1.x-=64.0/2;
	p2.x-=64.0/2;
	p1.z-=64.0/2;
	p2.z-=64.0/2;
	
//	float1 scl=15.0;
	float1 scl=45.0;
	p1*=scl;
	p2*=scl;
	
//	line(x3+p1.x,y3+p1.y,  x3+p2.x,y3+p2.y,color);
	line(x3+p1.x,y3+p1.z,  x3+p2.x,y3+p2.z,color);
}
void arrow3d(vec3 p1,vec3 p2,int color,float1 s=2.0)
{
	if(length(p2-p1)<delta_e) return;
	if(length(p2-p1)>1000.0) return;
	
	vec3 xaxis;
	vec3 yaxis;
	vec3 zaxis;
	vec3 base;
//coordinate axis of arrow	
	base=p1;
	yaxis=p2-p1;
//	float1 l=length(yaxis)-15.0;
	yaxis=normalize(yaxis);
	xaxis=normalize(cross(yaxis,vec3(yaxis.z,yaxis.x,yaxis.y)));
	zaxis=normalize(cross(xaxis,yaxis));


//arrow
	int u=12;//5
	vec3 p4;	
	line3d(p1,p2,color);
	line3d(p2,p2+(p1-p2)*0.5,color&0xff7f7f7f);//dark
//arrow head
//	float1 s=2.1;//2 0.1
if(needarrowhead)
	for(int i=0;i<u;i++)	
	{
		float1 aa=(float1)i*M_PI*2/u;
		p4=p2
			-xaxis*s*cos(aa)
			-zaxis*s*sin(aa)
			-yaxis*s*4.0;	
		line3d(p2,p4,color);
	}
}

void arrow3d(vec3 p1,vec3 p2,vec3 color)
{
	int col=0;
	col=(int)((color.x*255.0)); col<<=8;
	col+=(int)((color.y*255.0)); col<<=8;
	col+=(int)((color.z*255.0)); 
	arrow3d(p1,p2,col);
}

vec3 sphere(float1 alpha,float1 beta)
{
	alpha*=iradian;
	beta*=iradian;
	
	return vec3(
		sin(alpha)*cos(beta),
		sin(beta),
		cos(alpha)*cos(beta)	);
}

axis sphereaxis(float1 alpha,float1 beta)
{
	axis ax;
	vec3 p1=sphere(alpha,beta);
	vec3 p2=sphere(alpha+1e-3,beta);
	vec3 p3=sphere(alpha,beta+1e-3);
	
	ax.x=normalize(p2-p1);
	ax.z=normalize(p3-p1);
	ax.y=cross(ax.x,ax.z);
	
	return ax;
}
void calc_camaxis()
{
	camz=normalize(look-eye);
	camx=cross(camz,up);
	camy=cross(camx,camz);
	camx=normalize(camx);
	camy=normalize(camy);
}

float1 cos2(float1 a) 
{
	a=cos(a);
	return a*a;
}

float1 sin2(float1 a) 
{
	a=sin(a);
	return a*a;
}
inline int color(int r,int g,int b) 
{
	return (((r)<<16)+((g)<<8)+(b));
}

int colorv2int(vec3 col)
{
	int r=(int)(col.x*255.0);
	int g=(int)(col.y*255.0);
	int b=(int)(col.z*255.0);
	
	return (r<<16)+(g<<8)+b;
}

void rotvx(vec3 *ax,float1 alpha)
{
	alpha*=iradian;
	float1 c1=cos(alpha);
	float1 s1=sin(alpha);
	
	float1 x=ax[0].z*c1 + ax[0].y*s1;
	float1 y=ax[0].y*c1 - ax[0].z*s1;
	
	ax[0].z=x;
	ax[0].y=y;
}
void rotvy(vec3 *ax,float1 alpha)
{
	alpha*=iradian;
	float1 c1=cos(alpha);
	float1 s1=sin(alpha);
	
	float1 x=ax[0].z*c1 + ax[0].x*s1;
	float1 y=ax[0].x*c1 - ax[0].z*s1;
	
	ax[0].z=x;
	ax[0].x=y;
}
void rotvz(vec3 *ax,float1 alpha)
{
	alpha*=iradian;
	float1 c1=cos(alpha);
	float1 s1=sin(alpha);
	
	float1 x=ax[0].x*c1 + ax[0].y*s1;
	float1 y=ax[0].y*c1 - ax[0].x*s1;
	
	ax[0].x=x;
	ax[0].y=y;
}



double q_getcpuclock()
{
#ifdef __LINUX__
        timespec tim;
        timespec timres;
        clockid_t cid=0;

        clock_gettime(CLOCK_MONOTONIC,&tim);    

        double ti=(double)tim.tv_sec + (double)tim.tv_nsec/1e9;
        return ti;
#else
//	return (dd)GetTickCount()*0.001f;
        LARGE_INTEGER   freq;
        LARGE_INTEGER   counter;

        QueryPerformanceFrequency( &freq );
        QueryPerformanceCounter( &counter );
        
        return ( static_cast< double >( counter.QuadPart ) / 
                        static_cast< double >( freq.QuadPart ) );// * 100.0f
#endif
}
double ti22=0;

double xtimer(const char *ss)
{
	double ti1=q_getcpuclock();
	double deltaSeconds=ti1-ti22;
	ti22=ti1;

	if(deltaSeconds==0) deltaSeconds=0.000001;
	if(ss) printf("%s : %f sec %f fps\n",ss,deltaSeconds,1.0/deltaSeconds);
	return deltaSeconds;
}
int getcolor(float1 n,float1 h)
{
	float1 t1=n/h;
	int o=(int)t1;
	t1-=(float)o;
	
	vec3 col=lerp(coltab[o+1],coltab[o],t1);
	int coli=(((int)(col.x*255.0))<<16) + (((int)(col.y*255.0))<<8) + (((int)(col.z*255.0))) ;
	return coli;
}



struct qmat
{
	union { 
		struct { 
 			float1 _11,_12,_13,_14;
			float1 _21,_22,_23,_24;
			float1 _31,_32,_33,_34;
			float1 _41,_42,_43,_44;
		};
//		float1 m[16];
		float1 m[4][4];
	};
};


inline float determinant(const qmat *m)
{
	//amat m;memcpy(&m,m1,sizeof(qmat));
//11 12 13 14
//21 22 23 24
//31 32 33 34
//41 42 43 44 

	return
		m->_14 * m->_23 * m->_32 * m->_41
		-m->_13 * m->_24 * m->_32 * m->_41
		-m->_14 * m->_22 * m->_33 * m->_41
		+m->_12 * m->_24 * m->_33 * m->_41+
		
   		m->_13 * m->_22 * m->_34 * m->_41
   		-m->_12 * m->_23 * m->_34 * m->_41
   		-m->_14 * m->_23 * m->_31 * m->_42
   		+m->_13 * m->_24 * m->_31 * m->_42+
   		
   		m->_14 * m->_21 * m->_33 * m->_42
   		-m->_11 * m->_24 * m->_33 * m->_42
   		-m->_13 * m->_21 * m->_34 * m->_42
   		+m->_11 * m->_23 * m->_34 * m->_42+
   		
   		m->_14 * m->_22 * m->_31 * m->_43
   		-m->_12 * m->_24 * m->_31 * m->_43
   		-m->_14 * m->_21 * m->_32 * m->_43
   		+m->_11 * m->_24 * m->_32 * m->_43+
   		
   		m->_12 * m->_21 * m->_34 * m->_43
   		-m->_11 * m->_22 * m->_34 * m->_43
   		-m->_13 * m->_22 * m->_31 * m->_44
   		+m->_12 * m->_23 * m->_31 * m->_44+
   		
   		m->_13 * m->_21 * m->_32 * m->_44
   		-m->_11 * m->_23 * m->_32 * m->_44
   		-m->_12 * m->_21 * m->_33 * m->_44
   		+m->_11 * m->_22 * m->_33 * m->_44;
}
void qmatInverse(qmat *m1, qmat *m)//perfect!
{
	qmat n;
//11 12 13 14
//21 22 23 24
//31 32 33 34
//41 42 43 44 
	n._11 = m->_23*m->_34*m->_42 - m->_24*m->_33*m->_42 + m->_24*m->_32*m->_43 - m->_22*m->_34*m->_43 - m->_23*m->_32*m->_44 + m->_22*m->_33*m->_44;
	n._12 = m->_14*m->_33*m->_42 - m->_13*m->_34*m->_42 - m->_14*m->_32*m->_43 + m->_12*m->_34*m->_43 + m->_13*m->_32*m->_44 - m->_12*m->_33*m->_44;
	n._13 = m->_13*m->_24*m->_42 - m->_14*m->_23*m->_42 + m->_14*m->_22*m->_43 - m->_12*m->_24*m->_43 - m->_13*m->_22*m->_44 + m->_12*m->_23*m->_44;
	n._14 = m->_14*m->_23*m->_32 - m->_13*m->_24*m->_32 - m->_14*m->_22*m->_33 + m->_12*m->_24*m->_33 + m->_13*m->_22*m->_34 - m->_12*m->_23*m->_34;
	n._21 = m->_24*m->_33*m->_41 - m->_23*m->_34*m->_41 - m->_24*m->_31*m->_43 + m->_21*m->_34*m->_43 + m->_23*m->_31*m->_44 - m->_21*m->_33*m->_44;
	n._22 = m->_13*m->_34*m->_41 - m->_14*m->_33*m->_41 + m->_14*m->_31*m->_43 - m->_11*m->_34*m->_43 - m->_13*m->_31*m->_44 + m->_11*m->_33*m->_44;
	n._23 = m->_14*m->_23*m->_41 - m->_13*m->_24*m->_41 - m->_14*m->_21*m->_43 + m->_11*m->_24*m->_43 + m->_13*m->_21*m->_44 - m->_11*m->_23*m->_44;
	n._24 = m->_13*m->_24*m->_31 - m->_14*m->_23*m->_31 + m->_14*m->_21*m->_33 - m->_11*m->_24*m->_33 - m->_13*m->_21*m->_34 + m->_11*m->_23*m->_34;
	n._31 = m->_22*m->_34*m->_41 - m->_24*m->_32*m->_41 + m->_24*m->_31*m->_42 - m->_21*m->_34*m->_42 - m->_22*m->_31*m->_44 + m->_21*m->_32*m->_44;
	n._32 = m->_14*m->_32*m->_41 - m->_12*m->_34*m->_41 - m->_14*m->_31*m->_42 + m->_11*m->_34*m->_42 + m->_12*m->_31*m->_44 - m->_11*m->_32*m->_44;
	n._33 = m->_12*m->_24*m->_41 - m->_14*m->_22*m->_41 + m->_14*m->_21*m->_42 - m->_11*m->_24*m->_42 - m->_12*m->_21*m->_44 + m->_11*m->_22*m->_44;
	n._34 = m->_14*m->_22*m->_31 - m->_12*m->_24*m->_31 - m->_14*m->_21*m->_32 + m->_11*m->_24*m->_32 + m->_12*m->_21*m->_34 - m->_11*m->_22*m->_34;
	n._41 = m->_23*m->_32*m->_41 - m->_22*m->_33*m->_41 - m->_23*m->_31*m->_42 + m->_21*m->_33*m->_42 + m->_22*m->_31*m->_43 - m->_21*m->_32*m->_43;
	n._42 = m->_12*m->_33*m->_41 - m->_13*m->_32*m->_41 + m->_13*m->_31*m->_42 - m->_11*m->_33*m->_42 - m->_12*m->_31*m->_43 + m->_11*m->_32*m->_43;
	n._43 = m->_13*m->_22*m->_41 - m->_12*m->_23*m->_41 - m->_13*m->_21*m->_42 + m->_11*m->_23*m->_42 + m->_12*m->_21*m->_43 - m->_11*m->_22*m->_43;
	n._44 = m->_12*m->_23*m->_31 - m->_13*m->_22*m->_31 + m->_13*m->_21*m->_32 - m->_11*m->_23*m->_32 - m->_12*m->_21*m->_33 + m->_11*m->_22*m->_33;

    float1 det=determinant(m);

	float1 scale=1.0;
    if(det!=0.0) scale/=det;
    
	n._11*=scale;
	n._12*=scale;
	n._13*=scale;
	n._14*=scale;
	n._21*=scale;
	n._22*=scale;
	n._23*=scale;
	n._24*=scale;
	n._31*=scale;
	n._32*=scale;
	n._33*=scale;
	n._34*=scale;
	n._41*=scale;
	n._42*=scale;
	n._43*=scale;
	n._44*=scale;
	
	memcpy(m1,&n,sizeof(qmat));
}

inline void qmatMultiply(qmat *m1, const qmat *m2, const qmat *m3)
{
	qmat mo;

	mo._11=m2->_11*m3->_11 + m2->_12*m3->_21 + m2->_13*m3->_31 + m2->_14*m3->_41 ;
	mo._12=m2->_11*m3->_12 + m2->_12*m3->_22 + m2->_13*m3->_32 + m2->_14*m3->_42 ;
	mo._13=m2->_11*m3->_13 + m2->_12*m3->_23 + m2->_13*m3->_33 + m2->_14*m3->_43 ;
	mo._14=m2->_11*m3->_14 + m2->_12*m3->_24 + m2->_13*m3->_34 + m2->_14*m3->_44 ;

	mo._21=m2->_21*m3->_11 + m2->_22*m3->_21 + m2->_23*m3->_31 + m2->_24*m3->_41 ;
	mo._22=m2->_21*m3->_12 + m2->_22*m3->_22 + m2->_23*m3->_32 + m2->_24*m3->_42 ;
	mo._23=m2->_21*m3->_13 + m2->_22*m3->_23 + m2->_23*m3->_33 + m2->_24*m3->_43 ;
	mo._24=m2->_21*m3->_14 + m2->_22*m3->_24 + m2->_23*m3->_34 + m2->_24*m3->_44 ;

	mo._31=m2->_31*m3->_11 + m2->_32*m3->_21 + m2->_33*m3->_31 + m2->_34*m3->_41 ;
	mo._32=m2->_31*m3->_12 + m2->_32*m3->_22 + m2->_33*m3->_32 + m2->_34*m3->_42 ;
	mo._33=m2->_31*m3->_13 + m2->_32*m3->_23 + m2->_33*m3->_33 + m2->_34*m3->_43 ;
	mo._34=m2->_31*m3->_14 + m2->_32*m3->_24 + m2->_33*m3->_34 + m2->_34*m3->_44 ;

	mo._41=m2->_41*m3->_11 + m2->_42*m3->_21 + m2->_43*m3->_31 + m2->_44*m3->_41 ;
	mo._42=m2->_41*m3->_12 + m2->_42*m3->_22 + m2->_43*m3->_32 + m2->_44*m3->_42 ;
	mo._43=m2->_41*m3->_13 + m2->_42*m3->_23 + m2->_43*m3->_33 + m2->_44*m3->_43 ;
	mo._44=m2->_41*m3->_14 + m2->_42*m3->_24 + m2->_43*m3->_34 + m2->_44*m3->_44 ;

//	qmatCopy(m1,&mo);
	memcpy(m1,&mo,sizeof(qmat));
}

inline void qmatTranspose(qmat *m1,const qmat *m2)
{
	qmat mo;
	
	mo._11=m2->_11;
	mo._12=m2->_21;
	mo._13=m2->_31;
	mo._14=m2->_41;

	mo._21=m2->_12;
	mo._22=m2->_22;
	mo._23=m2->_32;
	mo._24=m2->_42;

	mo._31=m2->_13;
	mo._32=m2->_23;
	mo._33=m2->_33;
	mo._34=m2->_43;

	mo._41=m2->_14;
	mo._42=m2->_24;
	mo._43=m2->_34;
	mo._44=m2->_44;

	memcpy(m1,&mo,sizeof(qmat));
}



