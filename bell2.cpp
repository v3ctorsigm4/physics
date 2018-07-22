
#define __USE_X11
//https://arxiv.org/pdf/quant-ph/0402001


#include "eng6.cpp"



double S2classic(double a,double b,double a2,double b2)
{
	double N=0.0,step=M_PI/180.0,
		Pp1=0,Pm1=0,
		Pp2=0,Pm2=0,
		Pp3=0,Pm3=0,
		Pp4=0,Pm4=0;
	
	for(double p=0.0;p<2.0*M_PI;p+=step)
	{
		Pp1+=cos2(p-a)*cos2(p-b)*step;
		Pm1+=sin2(p-a)*cos2(p-b)*step;
	
		Pp2+=cos2(p-a)*cos2(p-b2)*step;
		Pm2+=sin2(p-a)*cos2(p-b2)*step;

		Pp3+=cos2(p-a2)*cos2(p-b)*step;
		Pm3+=sin2(p-a2)*cos2(p-b)*step;

		Pp4+=cos2(p-a2)*cos2(p-b2)*step;
		Pm4+=sin2(p-a2)*cos2(p-b2)*step;
	}
	Pp1/=(2.0*M_PI);	Pm1/=(2.0*M_PI);
	Pp2/=(2.0*M_PI);	Pm2/=(2.0*M_PI);
	Pp3/=(2.0*M_PI);	Pm3/=(2.0*M_PI);
	Pp4/=(2.0*M_PI);	Pm4/=(2.0*M_PI);
	
	N =(Pp1-Pm1)/(Pp1+Pm1);
	N-=(Pp2-Pm2)/(Pp2+Pm2);
	N+=(Pp3-Pm3)/(Pp3+Pm3);
	N+=(Pp4-Pm4)/(Pp4+Pm4);
	
	return N;
}
double S2retrocausal(double a,double b,double a2,double b2)
{
	double N=0.0,step=M_PI/180.0,
		Pp1=0,Pm1=0,
		Pp2=0,Pm2=0,
		Pp3=0,Pm3=0,
		Pp4=0,Pm4=0;
	
	for(double p=0.0;p<2.0*M_PI;p+=step)
	{
		for(int i=0;i<4;i++)
		{
			double w;

			if(i==0) w=a;
			if(i==1) w=b;
			if(i==2) w=a+M_PI/2;
			if(i==3) w=b+M_PI/2;
			Pp1+=cos2(w-a)*cos2(w-b)*step/4.0;
			Pm1+=sin2(w-a)*cos2(w-b)*step/4.0;
		
			if(i==0) w=a;
			if(i==1) w=b2;
			if(i==2) w=a+M_PI/2;
			if(i==3) w=b2+M_PI/2;
			Pp2+=cos2(w-a)*cos2(w-b2)*step/4.0;
			Pm2+=sin2(w-a)*cos2(w-b2)*step/4.0;
	
			if(i==0) w=a2;
			if(i==1) w=b;
			if(i==2) w=a2+M_PI/2;
			if(i==3) w=b+M_PI/2;
			Pp3+=cos2(w-a2)*cos2(w-b)*step/4.0;
			Pm3+=sin2(w-a2)*cos2(w-b)*step/4.0;
	
			if(i==0) w=a2;
			if(i==1) w=b2;
			if(i==2) w=a2+M_PI/2;
			if(i==3) w=b2+M_PI/2;
			Pp4+=cos2(w-a2)*cos2(w-b2)*step/4.0;
			Pm4+=sin2(w-a2)*cos2(w-b2)*step/4.0;
		}
	}
	Pp1/=(2.0*M_PI);	Pm1/=(2.0*M_PI);
	Pp2/=(2.0*M_PI);	Pm2/=(2.0*M_PI);
	Pp3/=(2.0*M_PI);	Pm3/=(2.0*M_PI);
	Pp4/=(2.0*M_PI);	Pm4/=(2.0*M_PI);
	
	N =(Pp1-Pm1)/(Pp1+Pm1);
	N-=(Pp2-Pm2)/(Pp2+Pm2);// - because a-b2 is the biggest difference
	N+=(Pp3-Pm3)/(Pp3+Pm3);
	N+=(Pp4-Pm4)/(Pp4+Pm4);
	
	return N;
}

double S2qm(double a,double b,double a2,double b2)
{
	double N=0.0,Pp=0,Pm=0;
	
	Pp=0.5*cos2(a-b);
	Pm=0.5*sin2(a-b);
	N=(Pp-Pm)/(Pp+Pm);

	Pp=0.5*cos2(a-b2);
	Pm=0.5*sin2(a-b2);
	N-=(Pp-Pm)/(Pp+Pm);

	Pp=0.5*cos2(a2-b);
	Pm=0.5*sin2(a2-b);
	N+=(Pp-Pm)/(Pp+Pm);

	Pp=0.5*cos2(a2-b2);
	Pm=0.5*sin2(a2-b2);
	N+=(Pp-Pm)/(Pp+Pm);

	return N;
}

double S2ww(double a,double b,double a2,double b2)
{
	double N=0.0,step=M_PI/180.0,
		Pp1=0,Pm1=0,
		Pp2=0,Pm2=0,
		Pp3=0,Pm3=0,
		Pp4=0,Pm4=0;
	
	double scl=2;//1.41;
	for(double p=0.0;p<2.0*M_PI;p+=step)
	{
		double p_wfield=frnd(M_PI*2.0);
//		double p_wfield=0.0;	if(frnd(1.0)>0.5) p_wfield=M_PI/2;
#if 0
		Pp1+=cos2(a-p_wfield)*cos2(b-p_wfield)*scl*step;
		Pm1+=sin2(a-p_wfield)*cos2(b-p_wfield)*scl*step;

		Pp2+=cos2(a-p_wfield)*cos2(b2-p_wfield)*scl*step;
		Pm2+=sin2(a-p_wfield)*cos2(b2-p_wfield)*scl*step;

		Pp3+=cos2(a2-p_wfield)*cos2(b-p_wfield)*scl*step;
		Pm3+=sin2(a2-p_wfield)*cos2(b-p_wfield)*scl*step;

		Pp4+=cos2(a2-p_wfield)*cos2(b2-p_wfield)*scl*step;
		Pm4+=sin2(a2-p_wfield)*cos2(b2-p_wfield)*scl*step;
#else
		Pp1+=cos2(a-p)*cos2(b-p)*cos2(a-p_wfield)*cos2(b-p_wfield)*scl*step;
		Pm1+=sin2(a-p)*cos2(b-p)*sin2(a-p_wfield)*cos2(b-p_wfield)*scl*step;

		Pp2+=cos2(a-p)*cos2(b2-p)*cos2(a-p_wfield)*cos2(b2-p_wfield)*scl*step;
		Pm2+=sin2(a-p)*cos2(b2-p)*sin2(a-p_wfield)*cos2(b2-p_wfield)*scl*step;

		Pp3+=cos2(a2-p)*cos2(b-p)*cos2(a2-p_wfield)*cos2(b-p_wfield)*scl*step;
		Pm3+=sin2(a2-p)*cos2(b-p)*sin2(a2-p_wfield)*cos2(b-p_wfield)*scl*step;

		Pp4+=cos2(a2-p)*cos2(b2-p)*cos2(a2-p_wfield)*cos2(b2-p_wfield)*scl*step;
		Pm4+=sin2(a2-p)*cos2(b2-p)*sin2(a2-p_wfield)*cos2(b2-p_wfield)*scl*step;
#endif
	}
	Pp1/=(2.0*M_PI);	Pm1/=(2.0*M_PI);
	Pp2/=(2.0*M_PI);	Pm2/=(2.0*M_PI);
	Pp3/=(2.0*M_PI);	Pm3/=(2.0*M_PI);
	Pp4/=(2.0*M_PI);	Pm4/=(2.0*M_PI);
	
	N =(Pp1-Pm1)/(Pp1+Pm1);
	N-=(Pp2-Pm2)/(Pp2+Pm2);
	N+=(Pp3-Pm3)/(Pp3+Pm3);
	N+=(Pp4-Pm4)/(Pp4+Pm4);
	
	return N;
}


void bell()
{
	init_system();
	
	line(0,100,180*4,100,0x555555);
	line(0,150,180*4,150,0x555555);
	line(0,200,180*4,200,0x999999);
	line(0,250,180*4,250,0x555555);
	line(0,300,180*4,300,0x555555);

	line(22*4,0,22*4,300,0x555555);//22.5 critical angle
	
    for(int angle=0;angle<180;angle+=2)
    {
	    double S=0.0,a,a2,b,b2;
	   
	    a=0.0*iradian;
	    b=angle*iradian;
	    a2=angle*2.0*iradian;
	    b2=angle*3.0*iradian;
	   

	    printf("%d:  ",angle);
	
	    S=S2classic(a,b,a2,b2 );    	printf("classic: %f ",S); 		pixel(angle*4,200-S*50.0,0x0000ff);//blue
	    S=S2retrocausal(a,b,a2,b2 );   	printf("retrocausal: %f  ",S); 	pixel(angle*4,200-S*50.0,0x00ff00);//green
    	S=S2qm(a,b,a2,b2 );    			printf("qm: %f ",S); 			pixel(angle*4+1,200-S*50.0,0xff0000);//red   
    	S=S2ww(a,b,a2,b2 );    			printf("WW: %f \n",S); 			pixel(angle*4+2,200-S*50.0,0xff00ff);//
	    qflush();
	}

    qflush();
    getchar();
}

int main()
{
	bell();
	return 0;
}


