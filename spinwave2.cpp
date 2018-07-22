
#define __USE_X11


#include "eng6.cpp"




void line2(vec3 v2,int color,int el)
{
	color = 255*(1-(v2.z)*0.01);
	if(v2.z<0) color *= 256*256;

	int x1 = 100+el*70;
	int y1 = 200;
	int x2 = x1 +v2.x;
	int y2 = 200+v2.z;
	
	XSetForeground(display5,gc,color);
	XDrawLine(display5, win, gc, x1,y1,x2,y2);
}
void point2(vec3 v1,int color,int el)
{
	int x = 100+v1.x+el*70;
	int y = 200+v1.z;

	XSetForeground(display5,gc,color);
	XDrawPoint(display5, win, gc, x,y);
}




void fnc()
{
	vec3 nor, mag, tmp;
	vec3 pos[20];
	vec3 vel[20];
	vec3 imp[20];
	float1 sz,dt = 0.01;
	mag=vec3(0,0,1);

	int count2 = 0,y,jj,kk,ll,N = 15;


	for( jj = 0;jj<N;jj++ )
	{
		float1 aa = 20*M_PI/180;

		pos[jj]=vec3( 30,0,0 );
		vel[jj]=vec3( 0,-cos(aa),sin(aa) );
		for(int ii=0;ii<3;ii++) vel[jj] *= 10.0;
	}
	for( jj = 0;jj<N/2;jj++ ) vel[jj].y*=-1.0;



	for( kk = 0;kk<3000000;kk++ )
	{
		for( jj = 0;jj<N;jj++ )
		{
			pos[jj] += vel[jj]*dt; //dot moving on sphere
			nor=normalize( pos[jj] );

			vel[jj] -= nor*dot( nor,vel[jj] );// cut normal component (covariant derivative)
			pos[jj]=normalize( pos[jj] )*30.0;//normalize

			vel[jj] += cross( vel[jj],mag )*0.1*dt;// Lorentz-like "force" (acceleration only)
			imp[jj]=normalize( cross( vel[jj],pos[jj] ) );// angular moment
		}
		for( jj = 0;jj<N;jj++ )
		for( ll = 0;ll<N;ll++ )
		if( jj != ll )
		{
			vel[jj] += cross( vel[jj],imp[ll] )*0.3*dt/fabs(10.0*sqr(jj-ll));//neigh acceleration
		}
		for( jj = 0;jj<N;jj++ )
			vel[jj] += cross( vel[jj],vec3(0,0,1) )*0.2*dt;//static field Z
			
		for( jj = 0;jj<N;jj++ )		vel[jj]=normalize( vel[jj] )* 10.0;


		count2++;
		if( count2>200 )
		{
			count2 = 0;
			qclear( );
		}


		for( jj = 0;jj<N;jj++ )
		{
			tmp = imp[jj]*100.0;
			line2( tmp,0xffff00,jj );
			point2( pos[jj],0x008f00,jj );

			if(imp[jj].z<0) XSetForeground( display5,gc,0xff0000 );
			else XSetForeground( display5,gc,0x0000ff );
			for(y = 0;y<8;y++)
			XDrawLine( display5, win, gc, 100+70*jj,300+y,100+70*jj+20,300+y );
		}
	}
}



int main()
{

	init_system();
	
	fnc();
	qflush();
	getchar();

	return 0;
}

