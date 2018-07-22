

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <X11/Xlib.h>
#include <assert.h>
#include <unistd.h>


Display *dpy;
Window w;
GC gc;

typedef double float1;
int ii;


void cross(float1 *v1,float1 *v2,float1 *v3)
{
	v1[0] = (v2[1]*v3[2]) - (v3[1]*v2[2]);
	v1[1] = (v2[2]*v3[0]) - (v3[2]*v2[0]);
	v1[2] = (v2[0]*v3[1]) - (v3[0]*v2[1]);
};
float1 dot(float1 *v1,float1 *v2) { return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];}

float1 length(float1 *v1) {return sqrt(dot(v1,v1));}

void normalize(float1 *v1,float1 *v2)
{
	float1 h= length(v2);
	for(int ii=0;ii<3;ii++) v1[ii] = v2[ii]/h;
}
float1 sqr(float1 v1) {return v1*v1;}
void def(float1 *v1,float1 x,float1 y,float1 z) {v1[0]=x; v1[1]=y; v1[2]=z;}




void line(float1 *v2,int color,int el)
{
	color = 255*(1-(v2[2])*0.01);
	if(v2[2]<0) color *= 256*256;

	int x1 = 100+el*70;
	int y1 = 200;
	int x2 = x1 +v2[0];
	int y2 = 200+v2[2];
	
	XSetForeground(dpy,gc,color);
	XDrawLine(dpy, w, gc, x1,y1,x2,y2);
}
void point(float1 *v1,int color,int el)
{
	int x = 100+v1[0]+el*70;
	int y = 200+v1[2];

	XSetForeground(dpy,gc,color);
	XDrawPoint(dpy, w, gc, x,y);
}




void fnc()
{
	float1 nor[3];
	float1 mag[3];
	float1 tmp[3];
	float1 pos[20][3];
	float1 vel[20][3];
	float1 imp[20][3];
	float1 sz,dt = 0.01;
	def(mag ,0,0,1);

	int count2 = 0,y,jj,kk,ll,N = 15;


	for( jj = 0;jj<N;jj++ )
	{
		float1 aa = 20*M_PI/180;

		def( pos[jj] , 30,0,0 );
		def( vel[jj] , 0,-cos(aa),sin(aa) );
		for(int ii=0;ii<3;ii++) vel[jj][ii] *= 10;
	}
	for( jj = 0;jj<N/2;jj++ ) vel[jj][1]*=-1.0;



	for( kk = 0;kk<3000000;kk++ )
	{
		for( jj = 0;jj<N;jj++ )
		{
			for(int ii=0;ii<3;ii++) pos[jj][ii] += vel[jj][ii]*dt;
			normalize( nor,pos[jj] );

			sz=dot( nor,vel[jj] );
			for(int ii=0;ii<3;ii++) vel[jj][ii] -= nor[ii]*sz;

			normalize( pos[jj],pos[jj] );
			for(int ii=0;ii<3;ii++) pos[jj][ii]*=30;


			cross( tmp,vel[jj],mag );
			for(int ii=0;ii<3;ii++) vel[jj][ii] += tmp[ii]*0.1*dt;

			cross( tmp,vel[jj],pos[jj] );
			normalize( imp[jj],tmp );
		}
		for( jj = 0;jj<N;jj++ )
		for( ll = 0;ll<N;ll++ )
		if( jj != ll )
		{
			cross( tmp,vel[jj],imp[ll] );
			for(int ii=0;ii<3;ii++) vel[jj][ii] += tmp[ii]*0.3*dt/fabs(10*sqr(jj-ll));
		}
		for( jj = 0;jj<N;jj++ )
		{
			normalize( vel[jj],vel[jj] );
			for(int ii=0;ii<3;ii++) vel[jj][ii] *= 10;
		}

		count2++;
		if( count2>200 )
		{
			count2 = 0;
			XClearWindow( dpy, w );
		}


		for( jj = 0;jj<N;jj++ )
		{
			for(int ii=0;ii<3;ii++) tmp[ii] = imp[jj][ii]*100;
			line( tmp,0xffff00,jj );
			point( pos[jj],0x008f00,jj );

			if(imp[jj][2]<0) XSetForeground( dpy,gc,0xff0000 );
			else XSetForeground( dpy,gc,0x0000ff );
			for(y = 0;y<8;y++)
			XDrawLine( dpy, w, gc, 100+70*jj,300+y,100+70*jj+20,300+y );
		}
	}
}



int main()
{
	dpy = XOpenDisplay((0));
	w = XCreateSimpleWindow(dpy, DefaultRootWindow(dpy), 0,0, 1200, 800, 0,0,0);

	XSelectInput(dpy, w, StructureNotifyMask);
	XMapWindow(dpy, w);

	gc = XCreateGC(dpy, w, 0, (0));
	XSetForeground(dpy,gc,0);

	for(;;) { XEvent e; XNextEvent(dpy, &e); if (e.type == MapNotify)break; }


	fnc();
	XFlush(dpy);
	getchar();

	return 0;
}

