//
//TETRAHEDRON of Stanrard Model  {weak hypercharge - weak isospin - electric charge}
//wikipedia: table of weak hypercharge 
//public domain
//g++ x.cpp -lX11 -O3

#define ENABLECONSOLEOUTPUT	
#define SHOWCUBE



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


#include <string>
#include <vector>
#include <map>

using namespace std;

Display *dpy;
Window win;
GC gc;


#define iradian (M_PI/180.0)
#define float1 float




float fsgn(float n)
{
	if(n<0.0) return -1.0;
	if(n>0.0) return 1.0;
	return 0.0;
}

struct vec3
{
float x,y,z;

vec3() {x=0;y=0;z=0;};
vec3(float x_,float y_,float z_) {x=x_;y=y_;z=z_;};

vec3 operator + (vec3 v) {vec3 v_;v_.x=x+v.x;v_.y=y+v.y;v_.z=z+v.z;return v_;};
vec3 operator - (vec3 v) {vec3 v_;v_.x=x-v.x;v_.y=y-v.y;v_.z=z-v.z;return v_;};
vec3 operator * (vec3 v) {vec3 v_;v_.x=x*v.x;v_.y=y*v.y;v_.z=z*v.z;return v_;};
vec3 operator / (vec3 v) {vec3 v_;v_.x=x/v.x;v_.y=y/v.y;v_.z=z/v.z;return v_;};

vec3 operator += (vec3 v) {x+=v.x;y+=v.y;z+=v.z;return *this;};
vec3 operator -= (vec3 v) {x-=v.x;y-=v.y;z-=v.z;return *this;};
vec3 operator *= (vec3 v) {x*=v.x;y*=v.y;z*=v.z;return *this;};
vec3 operator /= (vec3 v) {x/=v.x;y/=v.y;z/=v.z;return *this;};

vec3 operator + (float s) {vec3 v_;v_.x=x+s;v_.y=y+s;v_.z=z+s;return v_;};
vec3 operator - (float s) {vec3 v_;v_.x=x-s;v_.y=y-s;v_.z=z-s;return v_;};
vec3 operator * (float s) {vec3 v_;v_.x=x*s;v_.y=y*s;v_.z=z*s;return v_;};
vec3 operator / (float s) {vec3 v_;v_.x=x/s;v_.y=y/s;v_.z=z/s;return v_;};

vec3 operator += (float s) {x+=s;y+=s;z+=s;return *this;};
vec3 operator -= (float s) {x-=s;y-=s;z-=s;return *this;};
vec3 operator *= (float s) {x*=s;y*=s;z*=s;return *this;};
vec3 operator /= (float s) {x/=s;y/=s;z/=s;return *this;};

vec3 operator - () {x=-x;y=-y;z=-z;return *this;};
};

float dot(vec3 v1,vec3 v2) { return (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);}
float length(vec3 v1) {	return sqrt(dot(v1,v1));}
vec3 normalize(vec3 v1) { 	return v1/length(v1);}
vec3 cross(vec3 va,vec3 vb)
{
vec3 v_;

v_.x=(va.y*vb.z)-(vb.y*va.z);
v_.y=(va.z*vb.x)-(vb.z*va.x);
v_.z=(va.x*vb.y)-(vb.x*va.y);

return v_;
};

void cvariant(vec3 &v1) 	{	v1=v1/dot(v1,v1);}; 


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

float sqr(float n) {return n*n;}



int wx=1100;
int wy=650;

XEvent event;               
int mbut=0,mx,my,dmx,dmy,key=0;
vec3 camx,camy,camz,eye,look,up;

int transf3d2(vec3 v1,int *xc,int *yc)
{
	vec3 w1=v1-eye;
	float zoom=850;//650
	
	float x=dot(w1,camx);
	float y=dot(w1,camy);
	float z=dot(w1,camz);

	xc[0]=wx/2+ x*zoom/z;//perspective projection
	yc[0]=wy/2- y*zoom/z;
	if(z<0.0) return 0;
xc[0]-=150;
	
	return 1;
}
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

void line3d(vec3 p1,vec3 p2,int color)
{
	int x1,y1,x2,y2;

	if(transf3d2(p1,&x1,&y1))
	if(transf3d2(p2,&x2,&y2))
		line(x1,y1,x2,y2,color);
}
int colorv2int(vec3 col)
{
	int r=(int)(col.x*255.0);
	int g=(int)(col.y*255.0);
	int b=(int)(col.z*255.0);
	
	return (r<<16)+(g<<8)+b;
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

#if 1
    const char *pattern="*15*";//font manager   12,14,18
    int maxnames=10;
    int actual_count_return=0;
	char **fontnames=XListFonts(dpy, pattern, maxnames, &actual_count_return);
	if(actual_count_return==0)fontnames=XListFonts(dpy, "*14*", maxnames, &actual_count_return);
	if(actual_count_return==0)fontnames=XListFonts(dpy, "*18*", maxnames, &actual_count_return);
	if(actual_count_return==0)fontnames=XListFonts(dpy, "*10*", maxnames, &actual_count_return);
	if(actual_count_return==0)fontnames=XListFonts(dpy, "*12*", maxnames, &actual_count_return);
	if(actual_count_return==0)fontnames=XListFonts(dpy, "*13*", maxnames, &actual_count_return);
	for(int i=0;i<actual_count_return;i++)	printf("%s %d \n",fontnames[i],actual_count_return);
		
	XFontStruct *fontinfo=0; 
	fontinfo = XLoadQueryFont(dpy,fontnames[0]); 
	XSetFont(dpy,gc,fontinfo->fid); 
#endif



	for(;;) { XEvent e; XNextEvent(dpy, &e); if (e.type == MapNotify) break; }
	XMapWindow(dpy, win);


	return 0;
}
void qsleep(float dt)
{
	timespec ti1,ti2;

	ti1.tv_sec=0;
	ti1.tv_nsec=(int)(1e9*dt);

	nanosleep(&ti1,&ti2);
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





vector<string> str;//font
vector<int> posx;
vector<int> posy;
vector<int>	col;

float scale=600.0;

void printstr(char *buf,int x,int y,int color)
{
	XSetForeground(dpy,gc,color);
	XDrawString(dpy,win,gc,x,y,buf,strlen(buf));
}
void print3d(const char *buf,vec3 pos,int cols)
{
	int x=0,y=0;
	transf3d2(pos*scale,&x,&y);

	x+=3;//to right
	
	str.push_back(buf);
	posx.push_back(x);
	posy.push_back(y);
	col.push_back(cols);
}
void drawstrings()
{
	for(int i=0;i<str.size();i++)
		printstr((char *)str[i].c_str(),posx[i],posy[i],col[i]);
	str.clear();
	posx.clear();
	posy.clear();
	col.clear();
}



//the projection matrix
vec3 Yw_axis;//weak hypercharge
vec3 T3_axis;//weak isospin
vec3 Q_axis;//electric charge

int cursory;


void arrow(vec3 p1,vec3 p2,vec3 color,const char *mask,const char *mask2=0)
{	
//charge value dumps
	vec3 axial=p2-p1;
	float Yw=dot(axial,Yw_axis);//the projection matrix
	float T3=dot(axial,T3_axis);
	float Q =dot(axial,Q_axis);
	if(mask) 
	{
		char buf[164]={0};
#ifdef ENABLECONSOLEOUTPUT		
		printf("%s %s \t Q=%0.3f \t T3=%0.3f \t Yw=%0.3f \n",mask,mask2,Q,T3,Yw);
#endif
		sprintf(buf,"%s %s  Q=%0.3f  T3=%0.3f  Yw=%0.3f ",mask,mask2,Q,T3,Yw);
		printstr(buf,700,cursory,0xff00ff00);cursory+=15;

		print3d(mask,(p1+p2)/2.0,colorv2int(color));
	}
	
	
//draw	
	p1*=scale;
	p2*=scale;
	
	vec3 xaxis;
	vec3 yaxis;
	vec3 zaxis;
	vec3 base;
//coordinate axis of arrow	
	base=p1;
	yaxis=p2-p1;
	float l=length(yaxis)-15.0;
	yaxis=normalize(yaxis);
	xaxis=normalize(cross(yaxis,vec3(1,2,3)));
	zaxis=normalize(cross(xaxis,yaxis));


//arrow
	vec3 p4;	
	line3d(p1,p2,colorv2int(color));
//arrow head
	for(int i=0;i<20;i++)	
	{
		float aa=(float)i*M_PI*2/20.0;
		p4=p2
			-xaxis*5.0*cos(aa)
			-zaxis*5.0*sin(aa)
			-yaxis*20.0;	
		line3d(p2,p4,colorv2int(color));
	}
}


vec3 sphere5(float alpha,float beta)
{
	alpha*=iradian;
	beta*=iradian;
	
	vec3 v1(
		-sin(alpha)*cos(beta),
		sin(beta),
		cos(alpha)*cos(beta)	);

	vec3 yaxis=normalize(up);
	vec3 xaxis=normalize(cross(yaxis,vec3(yaxis.z,yaxis.x,yaxis.y)));
	vec3 zaxis=normalize(cross(xaxis,yaxis));
	
	vec3 v2=xaxis*v1.x + yaxis*v1.y + zaxis*v1.z;
	return v2;
}


void drawscene()
{
	vec3 p0,p1,p2,p3,p4,p5,p6,p62,p63		,a1,a2,c1,c2,c3,c4,cth1,cth2,cth3,cth4;
	
	printf("................................ \n");

//2 input parameters
	float edge=sqrt(3)/2; //length of real spin
	double r=0.5;//length of spin projection

	double edge2=cos(30.0*M_PI/180.0);//sqrt(3)/2   length of real spin
	double high=sqrt(edge2*edge2-r*r);
//	double alpha=acos(0.5/(sqrt(3.0)/2.0))/iradian;

	
	//cube base generator
	float cube_edge=sqrt(edge*edge/2);
	p1=vec3(-cube_edge/2.0,-cube_edge/2.0,-cube_edge/2.0);
	p2=vec3(-cube_edge/2.0, cube_edge/2.0, cube_edge/2.0);
	p3=vec3( cube_edge/2.0, cube_edge/2.0,-cube_edge/2.0);
	p4=vec3( cube_edge/2.0,-cube_edge/2.0, cube_edge/2.0);
	printf("%e \n",cube_edge);
	
	p0=(p1+p2+p3)/3.0;//face center
	p5=(p1+p2+p3+p4)/4.0;//volume center
	p6=(p2+p3+p4)/3.0;//face center
	p62=(p1+p2+p4)/3.0;//face center2
	p63=(p1+p3+p4)/3.0;//face center3




	cursory=10;


#ifdef SHOWCUBE
//cube +4 corner         vcenter + normalize(facecenter-vcenter)*length(p4-p5)
	c1=p5+ normalize(p0-p5)*length(p4-p5);
	c2=p5+ normalize(p6-p5)*length(p4-p5);
	c3=p5+ normalize(p62-p5)*length(p4-p5);
	c4=p5+ normalize(p63-p5)*length(p4-p5);
	// 1/3 of cube lines
	cth1=p3 + (c4-p3)/3.0;
	cth2=c2 + (p4-c2)/3.0;
	cth3=p2 + (c3-p2)/3.0;
	cth4=c1 + (p1-c1)/3.0;
#endif




//axis of charges
	vec3 centerofvL=(p3+p4)/2.0;
	vec3 centerofH=(p1+p2)/2.0;
	Q_axis=(centerofvL-centerofH);// Q -> orthogonal to vL   and H
	Yw_axis=(p4-p1)/2.0;// Yw  ->  direction of eR
	T3_axis=(p3-p2);// T3 -> direction of W
	
  
//	vec3 nor=normalize(cross(Yw_axis,T3_axis));	Yw_axis+=nor*0.1;	T3_axis+=nor*0.1;//nem jo, asszimetrikus lesz PL: eL vL

	cvariant(Yw_axis);//contravariant base vectors
	cvariant(T3_axis);
	cvariant(Q_axis);

printf("q %e \n",atan(Q_axis.z/Q_axis.x)/iradian);
printf("q %e \n",atan(Q_axis.x/Q_axis.y)/iradian);
printf("q %e \n",atan(Q_axis.y/Q_axis.z)/iradian);



//background
	int y=1;
	for(int x=-5;x<6;x++)	arrow(vec3(x,y,-6)/10.0,vec3(x,y,6)/10.0,vec3(0.1,0.1,0.1),0);
	for(int x=-5;x<6;x++)	arrow(vec3(-6,y,x)/10.0,vec3(6,y,x)/10.0,vec3(0.1,0.1,0.1),0);
	
//draw axis of charges
	float sc=0.4;
	arrow(p0,p0+Yw_axis*sc,vec3(0.9,0.9,0.9),0);
	arrow(p0,p0+T3_axis*sc,vec3(0.9,0.9,0.9),0);
	arrow(p0,p0+Q_axis*sc,vec3(0.9,0.9,0.9),0);
	print3d("Yw weak hypercharge",p0+Yw_axis*sc,0xffffff);
	print3d("T3 weak isospin",p0+T3_axis*sc,0xffffff);
	print3d("Q charge",p0+Q_axis*sc,0xffffff);




//draw particles	
//left handed  blue
	arrow(p4,p2,vec3(0.0,0.0,0.9),"eL","{4-2}");
	arrow(p4,p3,vec3(0.5,0.0,0.7),"vL","{4-3}");	
#if 0	
	a1=(p4+p3)/2.0;//vL
	a2=(p1+p2)/2.0;//H
	arrow(a1,a2,vec3(0.0,0.0,0.9),"eL3","{34-12}");
#endif
	arrow(p0,p3,vec3(0.0,0.0,0.5),"uL","{0-3}");
	arrow(p2,p6,vec3(0.0,0.0,0.5),"uL2","{2-6}");	
	arrow(p0,p2,vec3(0.0,0.0,0.5),"dL","{0-2}");	
	arrow(p3,p6,vec3(0.0,0.0,0.5),"dL2","{3-6}");	//alternative quarks

//right handed  red
	arrow(p4,p1,vec3(0.9,0.0,0.0),"eR","{4-1}");
	arrow(p0,p4,vec3(0.5,0.0,0.0),"uR","{0-4}");
	arrow(p1,p6,vec3(0.5,0.0,0.0),"uR2","{1-6}");

	arrow(p0,p1,vec3(0.5,0.0,0.0),"dR","{0-1}");
	arrow(p4,p6,vec3(0.5,0.0,0.0),"dR2","{4-6}");	


	a1=(p4+p1)/2.0;//eR
	a2=(p2+p3)/2.0;//W
	vec3 Z5=a2-a1;
	arrow(a1,a2,vec3(0.0,0.7,0.7),"Z","{41-23}");//and gamma
	arrow(p2,p3,vec3(0.9,0.9,0.0),"W","{2-3}");	

	arrow(p1,p2,vec3(0.1,0.6,0.0),"H","{1-2}");	
	arrow(p1,p3,vec3(0.3,0.5,0.0),"HQ","{1-3}");	

	a1=(p3+p1)/2.0;//1-H
	a2=(p4+p2)/2.0;//eL
//	arrow(a1,a2,vec3(0.0,0.9,0.0),"H2","{31-42}");//?


/*
	arrow(p1,p5,vec3(0.0,0.9,0.3),"W0","{1-5}");	
	arrow(p2,p5,vec3(0.0,0.9,0.3),"W1","{2-5}");	
	arrow(p3,p5,vec3(0.0,0.9,0.3),"W2","{3-5}");	
	arrow(p4,p5,vec3(0.0,0.9,0.3),"B0","{4-5}");	
*/
	float lenB=0.866025;
lenB=sqrt(cube_edge*cube_edge*3)/2; // /2
printf("lenB %e \n",lenB);
	vec3 pw3,pw1,pw2,pb0;
	pw3=normalize(p5-p1)*lenB;	arrow(p1,p1+pw3,vec3(0.0,0.9,0.3),"W3","{1-5}");	
	pw1=normalize(p5-p2)*lenB;	arrow(p2,p2+pw1,vec3(0.0,0.9,0.3),"W1","{2-5}");	
	pw2=normalize(p5-p3)*lenB;	arrow(p3,p3+pw2,vec3(0.0,0.9,0.3),"W2","{3-5}");	
	pb0=normalize(p5-p4)*lenB;	arrow(p4,p4+pb0,vec3(0.0,0.9,0.3),"B0","{4-5}");	




#ifdef SHOWCUBE
//cube
	vec3 cubecol=vec3(0.4,0.3,0.1)*0.5;
	arrow(p1,c3,cubecol,0,0);	
	arrow(c3,p4,cubecol,0,0);	
	arrow(p4,c4,cubecol,0,0);	
	arrow(c4,p1,cubecol,0,0);	

	arrow(c4,p3,cubecol,0,0);	
	arrow(p3,c1,cubecol,0,0);	
	arrow(c1,p1,cubecol,0,0);	

	arrow(c1,p2,cubecol,0,0);	
	arrow(p2,c2,cubecol,0,0);	
	arrow(p2,c2,cubecol,0,0);	

	arrow(p2,c3,cubecol,0,0);	
	arrow(c2,p3,cubecol,0,0);	
	arrow(c2,p4,cubecol,0,0);	

// 1/3 cube lines
	cubecol=vec3(0.1,0.4,0.3)*0.5;
	arrow(cth1,cth2,cubecol,0,0);	
	arrow(cth2,cth3,cubecol,0,0);	
	arrow(cth3,cth4,cubecol,0,0);	
	arrow(cth4,cth1,cubecol,0,0);	
#endif



//vertex numbers
	print3d("press 1-2-3-4",vec3(-0.5,0.0,-0.5),0xffffff);

	print3d("0",p0,0xffffff);
	print3d("1",p1,0xffffff);
	print3d("2",p2,0xffffff);
	print3d("3",p3,0xffffff);
	print3d("4",p4,0xffffff);
//	print3d("5",p5,0xffffff);	
	print3d("6",p6,0xffffff);	
	//print3d("8",p8,0xffffff);


	drawstrings();
	XFlush(dpy);//		qsleep(10);
}





int main()
{
	init_system();
	
	int first=1;
	float alpha=30,beta=30;
	up=normalize(vec3(1,0,0));//Q

	while(1)
	{
		if(keyboard() || first||key)
		{
			first=0;
			XClearWindow(dpy,win);
//camera		
		eye=sphere5(alpha,beta)*1000.0;
		look=vec3(0,0,0);
		eye+=look;

		alpha+=dmx*0.6;
		beta+=dmy*0.6;
//up vector		
		if(key=='1') up=normalize(vec3(1,0,0));//Q   
		if(key=='2') up=normalize(vec3(1,0,1));//Yw   
		if(key=='3') up=normalize(vec3(1,0,-1));//T3   
		if(key=='4') up=normalize(vec3(0,1,0));//norm  


		camz=normalize(look-eye);
		camx=cross(camz,up);
		camy=cross(camx,camz);
		camx=normalize(camx);
		camy=normalize(camy);
	

		drawscene();
		qsleep(100000.0);
		key=0;
		}
	}
	getchar();

	return 0;
}

#if 0	

eL {4-2} 	 Q=-1.000 	 T3=-0.500 	 Yw=-1.000 
vL {4-3} 	 Q=0.000 	 T3=0.500 	 Yw=-1.000 

uL {0-3} 	 Q=0.667 	 T3=0.500 	 Yw=0.333 
uL2 {2-6} 	 Q=0.667 	 T3=0.500 	 Yw=0.333 
dL {0-2} 	 Q=-0.333 	 T3=-0.500 	 Yw=0.333 
dL2 {3-6} 	 Q=-0.333 	 T3=-0.500 	 Yw=0.333 

eR {4-1} 	 Q=-1.000 	 T3=0.000 	 Yw=-2.000 
uR {0-4} 	 Q=0.667 	 T3=0.000 	 Yw=1.333 
uR2 {1-6} 	 Q=0.667 	 T3=0.000 	 Yw=1.333 
dR {0-1} 	 Q=-0.333 	 T3=0.000 	 Yw=-0.667 
dR2 {4-6} 	 Q=-0.333 	 T3=0.000 	 Yw=-0.667 

Z {41-23} 	 Q=0.000 	 T3=0.000 	 Yw=0.000 
W {2-3} 	 Q=1.000 	 T3=1.000 	 Yw=0.000 
H {1-2} 	 Q=0.000 	 T3=-0.500 	 Yw=1.000 
HQ {1-3} 	 Q=1.000 	 T3=0.500 	 Yw=1.000 

#endif


