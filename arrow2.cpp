
#define __USE_X11


#include "eng6.cpp"


#define _DU_DIRECTION


vec3 sphere5(float1 u,float1 v,float1 r)
{
	return vec3(
		r*cos(u)*cos(v),
		r*sin(v),
		r*sin(u)*cos(v));
}
void drawscene()
{
	float1 r=300.0;
	float1 du=M_PI/30.0;
	float1 dv=M_PI/30.0;

	for(int i=0;i<60;i++)
	for(int j=0;j<30;j++)
	{
		float1 u=i*du;
		float1 v=-M_PI/2+j*dv;
		vec3 p1=sphere5(u,v,r);
		vec3 p2=sphere5(u+du,v,r);
		vec3 p3=sphere5(u   ,v+dv,r);
		vec3 normal=cross(p2-p1,p3-p1);

		if(dot(p1-eye,normal)>0.0)
		{
			arrow3d(p1,p2,vec3(1.0,0.4,0.0)*0.4);	
			arrow3d(p1,p3,vec3(0.0,0.4,1.0)*0.4);	
		}
	
	}

	vec3 P0(-1,-1,-1);
	
int i=8,j=20;
#ifdef _DU_DIRECTION
	for(i=5;i<10;i++)// du
#else		
	for(j=20;j<25;j++)// dv
#endif		
	{
		float1 u=i*du;
		float1 v=-M_PI/2+j*dv;
		vec3 P1 =sphere5(u,v,r);
		vec3 P2 =sphere5(u+du,v,r);
		vec3 P3 =sphere5(u   ,v+dv,r);
		vec3 Pu=(P2-P1)/du;// u first derivative
		vec3 Pv=(P3-P1)/dv;// v 
		vec3 normal=normalize(cross(Pv,Pu));//surface normal
		
		if(P0.z==-1) P0=P1;//base point for second derivatives
		

		u+=du;//P1 moved into u direction   P1b
		vec3 P1b =sphere5(u,v,r);
		vec3 P2b =sphere5(u+du,v,r);
		vec3 P3b =sphere5(u   ,v+dv,r);
		vec3 Pub=(P2b-P1b)/du;
		vec3 Pvb=(P3b-P1b)/dv;

		u-=du; //restore
		v+=dv;//P1 moved into v direction     P1c
		vec3 P1c =sphere5(u,v,r);
		vec3 P2c =sphere5(u+du,v,r);
		vec3 P3c =sphere5(u   ,v+dv,r);
		vec3 Puc=(P2c-P1c)/du;
		vec3 Pvc=(P3c-P1c)/dv;
		
		vec3 Puu=(Pub-Pu)/du;// second derivatives
		vec3 Puv=(Puc-Pu)/dv;
		vec3 Pvu=(Pvb-Pv)/du;
		vec3 Pvv=(Pvc-Pv)/dv;

//first derivatives
#if 0
		arrow3d(P1,P1+Pu,vec3(1.0,0.0,0.0));	//at original place
		arrow3d(P1,P1+Pv,vec3(0.0,1.0,0.0));	
		arrow3d(P1,P1+normal*100.0,vec3(0.0,0.0,1.0));	
#else
		arrow3d(P0,P0+Pu,vec3(1.0,0.0,0.0));	//all from "zero" point
		arrow3d(P0,P0+Pv,vec3(0.0,1.0,0.0));	
		arrow3d(P0,P0+normal*100.0,vec3(0.0,0.0,1.0));	
#endif		
		
		float1 E=dot(Pu,Pu);//g(uu)         element of metric tensor
		float1 F=dot(Pu,Pv);//g(uv)=(gvu)
		float1 G=dot(Pv,Pv);//g(vv)
		
//contravariant of tangents  -> because the tangent is contravariant already, this is the covarian tangent xD		
		vec3 Pu_= Pu/E;// == Pu/dot(Pu,Pu) == normalize(Pu)/length(Pu)   because E==length(Pu)^2
		vec3 Pv_= Pv/G;
		
		float1 Christoffel_u_uu=dot(Puu,Pu_);
		float1 Christoffel_u_uv=dot(Puv,Pu_);
		float1 Christoffel_u_vu=dot(Pvu,Pu_);
		float1 Christoffel_u_vv=dot(Pvv,Pu_);
		
		float1 Christoffel_v_uu=dot(Puu,Pv_);
		float1 Christoffel_v_uv=dot(Puv,Pv_);
		float1 Christoffel_v_vu=dot(Pvu,Pv_);
		float1 Christoffel_v_vv=dot(Pvv,Pv_);

		float1 normalcomp_uu=dot(Puu,normal);//L
		float1 normalcomp_uv=dot(Puv,normal);//M
		float1 normalcomp_vu=dot(Pvu,normal);
		float1 normalcomp_vv=dot(Pvv,normal);//N


//second derivatives
#ifdef _DU_DIRECTION
		arrow3d(P0+Pu,P0+Pu+Puu*du,vec3(1.0,0.0,1.0));
		arrow3d(P0+Pv,P0+Pv+Pvu*du,vec3(0.0,1.0,1.0));
//Gauss–Codazzi solution
		vec3 PuuGC=Pu*Christoffel_u_uu + Pv*Christoffel_v_uu + normal*normalcomp_uu; // == Puu
		arrow3d(P0+Pu,P0+Pu + PuuGC*du,vec3(1.0,1.0,1.0));
#else
		arrow3d(P0+Pu,P0+Pu+Puv*dv,vec3(1.0,0.0,1.0));
		arrow3d(P0+Pv,P0+Pv+Pvv*dv,vec3(0.0,1.0,1.0));
//Gauss–Codazzi solution
		vec3 PuvGC=Pu*Christoffel_u_uv + Pv*Christoffel_v_uv + normal*normalcomp_uv; // == Puv
		arrow3d(P0+Pu,P0+Pu + PuvGC*du,vec3(1.0,1.0,1.0));
#endif


	}
}



int main()
{
	init_system();
	
	eye=vec3(400,600,400);
	look=vec3(0,0,0);
	up=vec3(0,1,0);
	calc_camaxis();

	qclear();
	drawscene();

	qflush();
	getchar();
	
	return 0;
}



