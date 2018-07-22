

//__________________________________________________________________

#define __USE_X11


#include "eng6.cpp"



//__________________________________________________________________
//__________________________________________________________________
//__________________________________________________________________
//__________________________________________________________________


//reference(400s):interf:200,  QWP(400s):80,   ERASER(800s):120/2=60       



#define ERASER_PRESENT 		1
#define QWP_PRESENT 		2
#define CLASSIC_THEORY 		4
#define TIMEREVERSE_THEORY 	8



//
//int setup=0 ;
//int setup=QWP_PRESENT ;
int setup=ERASER_PRESENT | QWP_PRESENT ;



int dslitepr()
{
	init_system();

	setup|= CLASSIC_THEORY;

	double scale=580.0;//visual scale (400s)
scale=400.0;
scale=540.0;
scale=680.0;				//80
	if(setup&ERASER_PRESENT) scale*=2.0; //(800s)


	for(int i=100;i<500;i+=20)      line(0,i,400,i,0x005500);
        line(0,500,400,500,0x00aa00);
        line(0,300,400,300,0x00aa00);
        line(0,100,400,100,0x00aa00);

    
    double eraser_angle=(-45.0*1.0)*iradian;//+-45
//eraser_angle=(0.0*1.0)*iradian;//+-45

    for(int ds_x=0;ds_x<200;ds_x++)//Ds position  -+4mm    
    {
        double P=0,Pds=0,Pdp=0;
        
        double step=M_PI*2.0/300.0;
        for(double p=0;p<M_PI*2.0;p+=step)
        {
            vec2 amp_dp,amp_ds,amp2;
            double photon_pol_a, photon_pol_b;


            {
            	double p2=p;//  almost good classic solution:   p2=p 
	            photon_pol_a=p2;	           	photon_pol_b=p2;//+M_PI/2;
//	            if(rand()%100<50)    {	photon_pol_a=p2+M_PI/2;photon_pol_b=p2;    }
	        }
     
            double dphase1=M_PI*frand();
//			double dphase2=dphase1;if(rand()%100>50) dphase2=M_PI*frand()*2.0;
			

            double ds_distance=1250.0-420.0;//mm  125-42 cm
            double dp_distance=980.0; //98 cm
            double wavelength1=702.2e-6;//mm      e-9m
            double k1=2.0*M_PI/wavelength1;//wave number
            double wavelength2=452.2e-6;
            double k2=2.0*M_PI/wavelength2;//wave number
            double hole_dist=0.2;//0.2
            double hole_wide=0.2;//0.2mm  = 200 micrometer wide   
            double ds_pos=4.0*(double)(ds_x-100)/100.0;//+-4mm     Ds position
            double dist3=dp_distance;
            int slit_wide=10;
			
            for(int w=0;w<slit_wide;w++)
            {
                double hole1x=hole_dist/2.0 + hole_wide*(double)w/slit_wide;//hole A
                double hole2x=-hole_dist/2.0 - hole_wide*(double)w/slit_wide;//hole B
                double dist1=sqrt(sqr(ds_pos - hole1x) + sqr(ds_distance));
                double dist2=sqrt(sqr(ds_pos - hole2x) + sqr(ds_distance));    

                add_amp(&amp_dp,photon_pol_a  ,1.0);//DP side

                if(setup&QWP_PRESENT)
                {
                   	add_quarterwaveplatel2c(&amp_ds,photon_pol_b, dist1*k1 +dphase1,-45.0*iradian,0.5);// DS side
                    add_quarterwaveplatel2c(&amp_ds,photon_pol_b, dist2*k2 +dphase1, 45.0*iradian,0.5);
                }   
                else
                {
                    add_amp(&amp_ds,photon_pol_b +dist1*k1 +dphase1,0.5);//2 slit  DS side
                    add_amp(&amp_ds,photon_pol_b +dist2*k2 +dphase1,0.5);
  //                  add_polarizer2(&amp_ds,45.0*iradian);
                }
            }
            double normalization_factor=1.0/(double)slit_wide;
            amp_ds=amp_ds*normalization_factor;
            amp_dp=amp_dp*normalization_factor;
            
            
            if(setup&ERASER_PRESENT)
                    add_polarizer2(&amp_dp,eraser_angle);// polarizer(eraser) before Dp  

            P+=probability(amp_dp,amp_dp)*probability(amp_ds,amp_ds)*step;
            Pds+=probability(amp_ds,amp_ds)*step;
            Pdp+=probability(amp_dp,amp_dp)*step;
        }
        double normalization_factor2=(M_PI*2.0);
        P/=normalization_factor2;       pixel(ds_x*2,501-P*scale,0xffff00);
        Pds/=normalization_factor2;     pixel(ds_x*2,500-Pds*scale,0xff0000);
        Pdp/=normalization_factor2;     pixel(ds_x*2,502-Pdp*scale,0x00ff);
    }

    //XFlush(dpy);    getchar();
        
    return 0;
}


int dslitepr2()
{
	init_system();

	setup ^= CLASSIC_THEORY;	setup |= TIMEREVERSE_THEORY;
//	setup|= CLASSIC_THEORY;

	double scale=580.0;//visual scale (400s)
scale=530.0;	
	if(setup&ERASER_PRESENT) scale*=2.0; //(800s)

	int dx=500;
	for(int i=100;i<500;i+=20)      line(dx,i,dx+400,i,0x005500);
        line(dx,500,dx+400,500,0x00aa00);
        line(dx,300,dx+400,300,0x00aa00);
        line(dx,100,dx+400,100,0x00aa00);

    
    double eraser_angle=(45.0*1.0)*iradian;//+-45
//eraser_angle=(-0.0*1.0)*iradian;//+-45

    for(int ds_x=0;ds_x<200;ds_x++)//Ds position  -+4mm    
    {
        double P=0,Pds=0,Pdp=0;
        
        double step=M_PI*2.0/300.0;
        for(double p=0;p<M_PI*2.0;p+=step)
        {
            vec2 amp_dp,amp_ds;
            double photon_pol_a, photon_pol_b;

/////            if(setup&CLASSIC_THEORY)
            {
            	double p2=0;//  almost good classic solution:   p2=p 
	            photon_pol_a=p2;	           	photon_pol_b=p2+M_PI/2;
//	            if(rand()%100<50)    {	photon_pol_a=p2+M_PI/2;photon_pol_b=p2;    }
	        }
            if(setup&TIMEREVERSE_THEORY)
			{
	            if(setup&ERASER_PRESENT)
	            {
	            	photon_pol_a=eraser_angle;//source modified by Dp polarizer (eraser)  RIGHT solution
		           	photon_pol_b=eraser_angle+M_PI/2;

		            if(rand()%100<50)    {	photon_pol_a=eraser_angle+M_PI/2;photon_pol_b=eraser_angle;    }
		        }
	        }
     

            double dphase1=M_PI*2*frand();
//			double dphase2=dphase1;if(rand()%100>50) dphase2=M_PI*frand()*2.0;
			
//something destroy the half of sample!! no matter what it is must to emulate
//if(rand()%100>40) dphase2=2*M_PI*frand();


            double ds_distance=1250.0-420.0;//mm  125-42 cm
            double dp_distance=980.0; //98 cm
            double wavelength=702.2e-6;//mm      e-9m
            double k1=2.0*M_PI/wavelength;//wave number
          double wavelength2=452.2e-6;
            double k2=2.0*M_PI/wavelength2;//wave number
              double hole_dist=0.2;//0.2
            double hole_wide=0.2;//0.2mm  = 200 micrometer wide   
            double ds_pos=4.0*(double)(ds_x-100)/100.0;//+-4mm     Ds position
            int slit_wide=10;
			
            for(int w=0;w<slit_wide;w++)
            {
                double hole1x=hole_dist/2.0 + hole_wide*(double)w/slit_wide;//hole A
                double hole2x=-hole_dist/2.0 - hole_wide*(double)w/slit_wide;//hole B
                double dist1=sqrt(sqr(ds_pos - hole1x) + sqr(ds_distance));
                double dist2=sqrt(sqr(ds_pos - hole2x) + sqr(ds_distance));    

                if(setup&QWP_PRESENT)
                {
	                add_amp(&amp_dp,photon_pol_a ,1.0);//DP side
                   	add_quarterwaveplatel2c(&amp_ds,photon_pol_b, dist1*k1 +dphase1,-45.0*iradian,0.5);//2 slit with QWP   DS side
                    add_quarterwaveplatel2c(&amp_ds,photon_pol_b, dist2*k2 +dphase1, 45.0*iradian,0.5);
//          amp_dp=amp_ds;rot(amp_dp);normalize(amp_dp);//infinite speed test  FALSE
          
                }   
                else
                {
	                add_amp(&amp_dp,photon_pol_a ,1.0);//DP side
                    add_amp(&amp_ds,photon_pol_b +dist1*k1 +dphase1,0.5);//2 slit  DS side
                    add_amp(&amp_ds,photon_pol_b +dist2*k2 +dphase1,0.5);
// virtual polarizer(eraser) after slits     test: polarization change AFTER QWP and slits(before Ds detector)
//            add_polarizer2(&amp_ds,eraser_angle);//FALSE
                }
            }
            double normalization_factor=1.0/(double)slit_wide;
            amp_ds=amp_ds*normalization_factor;
            amp_dp=amp_dp*normalization_factor;
            
            
            if(setup&ERASER_PRESENT)
                    add_polarizer2(&amp_dp,eraser_angle);// polarizer(eraser) before Dp  

            P+=probability(amp_dp,amp_dp)*probability(amp_ds,amp_ds)*step;
            Pds+=probability(amp_ds,amp_ds)*step;
            Pdp+=probability(amp_dp,amp_dp)*step;
        }
        double normalization_factor2=(M_PI*2.0);
        P/=normalization_factor2;       pixel(dx+ds_x*2,501-P*scale,0xffff00);
        Pds/=normalization_factor2;     pixel(dx+ds_x*2,500-Pds*scale,0xff0000);
        Pdp/=normalization_factor2;     pixel(dx+ds_x*2,502-Pdp*scale,0x00ff);
    }

    qflush();
    getchar();
        
    return 0;
}



int qwp()
{
	init_system();


	double r=50.0;
	for(double i=0;i<M_PI*2*3;i+=0.001)
	{
		vec2 amp; 							
		add_quarterwaveplatel2c(&amp,30.0*iradian, i ,30.0*iradian,1.0);
       	int x=150+amp.x*r;
       	int y=250-amp.y*r;
       	pixel(x,y,0xffff00);

      	r+=0.005;
	}
	r=50.0;
	for(double i=0;i<M_PI*2*3;i+=0.001)
	{
		vec2 amp; 							
		add_quarterwaveplatel2c(&amp,30.0*iradian, i ,(30.0+45.0)*iradian,1.0);
       	int x=350+amp.x*r;
       	int y=250-amp.y*r;
       	pixel(x,y,0xffff00);
      	r+=0.005;
	}
	r=50.0;
	for(double i=0;i<M_PI*2*3;i+=0.001)
	{
		vec2 amp; 							
		add_quarterwaveplatel2c(&amp,30.0*iradian, i ,(30.0-45.0)*iradian,1.0);
       	int x=550+amp.x*r;
       	int y=250-amp.y*r;
       	pixel(x,y,0xffff00);

      	r+=0.005;
	}
	r=50.0;
	for(double i=0;i<M_PI*2*3;i+=0.001)
	{
		vec2 amp; 							
		add_quarterwaveplatel2c(&amp,30.0*iradian, i ,(30.0-21.0)*iradian,1.0);
       	int x=750+amp.x*r;
       	int y=250-amp.y*r;
       	pixel(x,y,0xffff00);

      	r+=0.005;
	}

    qflush();
    getchar();
        
    return 0;
}


/*double cos2(double a) 
{
	a=cos(a);
	return a*a;
}*/


int epr()
{
	init_system();

	line(0,200, 400,200,0x555555);
	line(0,400, 400,400,0x555555);
	 
	double b=22.5*iradian;
    for(int x=0;x<360;x++)
    {
        double P=0,P2=0,P3=0,a=iradian*x,step=M_PI*2.0/1000.0;

	    for(double p=0;p<M_PI*2.0;p+=step)
    	{
	        P+=cos2(p-a)*cos2(p-b)*step;
        }
        P/=(M_PI*2.0);
        pixel(x,400-(int)(P*200.0),0x0000ff);//blue classic
        

#if 1
		P=0;
	    for(double p=0;p<M_PI*2.0;p+=step)
    	{
    		vec2 amp1,amp2;
    		double t2=frand()*M_PI*2.0;//(double)(rand()%100000)/100.0;
t2*=(5.0/5.0);
//t2+=10000000.0;

    		add_amp(&amp1,-t2/(5.0),0.5);
    		add_amp(&amp1, t2/(4.9),0.5);
			add_polarizer2(&amp1,a);
			double P1=probability(amp1,amp1);
			
    		add_amp(&amp2,-t2/(5.0),0.5);
    		add_amp(&amp2, t2/(4.9),0.5);
			add_polarizer2(&amp2,b);
			double P2=probability(amp2,amp2);

	        P+=P2*P1*step;
        }
        P/=(M_PI*2.0);
        pixel(x,400-(int)(P*200.0),0x00ff00);//blue classic
        
#else
		P=0;
	    for(double p=0;p<M_PI*2.0;p+=step)
    	{
	        P+=cos2(a-a)*cos2(a-b)*step*0.25;
	        P+=cos2(b-a)*cos2(b-b)*step*0.25;
	        P+=cos2(a+90.0*iradian -a)*cos2(a+90.0*iradian -b)*step*0.25;
	        P+=cos2(b+90.0*iradian -a)*cos2(b+90.0*iradian -b)*step*0.25;
        }
        P/=(M_PI*2.0);
        pixel(x+1,400-(int)(P*200.0),0x00ff00);//green retrocausal
#endif

        P=cos2(a-b)*0.5;
        pixel(x,400-(int)(P*200.0),0xff0000);// red QM
    }


    qflush();
    getchar();
        
    return 0;
}
//__________________________________________________________________
#if 1
void elec()
{
    float c=3e8,h=6.626e-34,m=9.1e-31,scale=1e-12,
	w,v,cf,fi1,gamma,l1,l4d,f,f1,f2,lm,l1a,l1b,fia,fib,fa;
    vec2 ewave_dir,lightcone_dir1,lightcone_dir2;
    float E1,E2,f3,f4;
    
    scale=h/(m*c)/50.0;
//  scale=h/(m*c)/20.0;
//  scale=h/(m*c)/10.0;

	init_system();


#if 1  
	double f1a,f1b,f2a,f2b;  
    v=0.23*c;
    cf=c*c/v;//QM phase velocity
    fi1=atan(v/c);//angle in spacetime
    gamma=1.0/sqrt(1.0-v*v/(c*c));
    l1=h/(m*gamma*v);//De Broglie wavelength
    E1=m*gamma*c*c;

#if 1
    f=m*c*c*gamma/h;//frequency
    f1a=f*((c+v)/(c));//Classic Doppler
    f1b=f*((c-v)/(c));
#else    
    l4d=l1*sin(fi1);//spacetime  wavelength
    lm=l4d/cos(fi1*2.0);    f=sqrt(c*c+v*v)/lm;//        f=f*(cf/(cf-v));
    f1a=f*(c/(c-v));//Classic Doppler
    f1b=f*(c/(c+v));//kisebb
#endif

    v=0.476*c;
    cf=c*c/v;//QM phase velocity
    fi1=atan(v/c);//angle in spacetime
    gamma=1.0/sqrt(1.0-v*v/(c*c));
    l1=h/(m*gamma*v);//De Broglie wavelength
    E2=m*gamma*c*c;

#if 1
    f=m*c*c*gamma/h;//frequency
    f2a=f*((c+v)/(c));//Classic Doppler
    f2b=f*((c-v)/(c));
#else    
    l4d=l1*sin(fi1);//spacetime  wavelength
    lm=l4d/cos(fi1*2.0);    f=sqrt(c*c+v*v)/lm;//        f=f*(cf/(cf-v));
    f2a=f*(c/(c-v));//Classic Doppler
    f2b=f*(c/(c+v));//kisebb
#endif

//a foton frekije egyszerre fugg a negy komponenstol /2 electron/    
//    f3=((f2a-f1a)+(f2b-f1b))/2;

//f1b*=sin(45.0*iradian);f2b*=sin(45.0*iradian);

    f3=(f2b-f1b)/1;
    f4=(E2-E1)/h;
printf("%e %e %e ?\n",f3,f4,f3/f4);
#endif


    
    v=0.23*c;
    cf=c*c/v;//QM phase velocity
    fi1=atan(v/c);//angle in spacetime
    gamma=1.0/sqrt(1.0-v*v/(c*c));
    l1=h/(m*gamma*v);//De Broglie wavelength
    
    l4d=l1*sin(fi1);//spacetime  wavelength
    ewave_dir.x=sin(fi1);// direction in spacetime
    ewave_dir.y=cos(fi1);
    f=m*c*c*gamma/h;//frequency
    fa=f;
//    f=m*c*c/h;//frequency
    w=f*2.0*M_PI;//angular frequency 


#if 1
    f1=fa*((c+v)/c);//Classic Doppler
    f2=fa*((c-v)/c);//kisebb
    f=(f1+f2)/2.0;
printf("fw %e \n",f);
#else
printf("fa %e \n",fa);
double f0=m*c*c/h;
#if 1
#else
f=f0*((cf-v)/(cf-v));printf("f1 %e \n",f);
f=f0*((cf+v)/(cf+v));printf("f2 %e \n",f);
f=f0*(cf/(cf-v));printf("fx %e \n",f);
f=f0*(cf/(cf+v));printf("fx %e \n",f);
f=f0*((cf-v)/cf);printf("fc %e \n",f);//f0=f
f=f0*((cf+v)/cf);printf("fx %e \n",f);
#endif
    lm=l4d/cos(fi1*2.0);    f=sqrt(c*c+v*v)/lm;//        f=f*(cf/(cf-v));
printf("fb %e \n",f);

    f1=f0*(sqrt(1.0+v/c) / sqrt(1.0-v/c));//OKKKK!!!!!!!!!!!!!!!1
    f2=f0*(sqrt(1.0-v/c) / sqrt(1.0+v/c));
//    f1=f*(c/(c-v));//Classic Doppler
//    f2=f*(c/(c+v));//kisebb
#endif
    f=(f1+f2)/2.0;
printf("fq %e \n",f);
printf("fa %e \n",fa);


    gamma=1.0/sqrt(1.0-v*v/(c*c));//timedilation
printf("gamma %e \n",gamma);
	gamma=c/sqrt(c*c-v*v);//lightclock
printf("gamma %e \n",gamma);


//f1 jobbra dolo ferde
    l1a=c/f1*sin(45.0*iradian);//a component on lightcone
    l1b=c/f2*sin(45.0*iradian);//b component on lightcone
    
 

    
    fia=45.0*iradian;
    lightcone_dir1.x=sin(fia);
    lightcone_dir1.y=cos(fia);
    fib=-45.0*iradian;
    lightcone_dir2.x=sin(fib);
    lightcone_dir2.y=cos(fib);

    
    
    for(int y=0;y<550;y++)
    for(int x=0;x<800;x++)
    {
	    int red=0,green=0;
		float amp1=0,phase,l,k;
		vec2 amp,screen;
        
		screen.x=scale*x;
		screen.y=scale*y;

#if 1
//spacetime wave      4d wavevector
		k=M_PI*2.0/l4d;
		phase=dot(screen,ewave_dir)*k;
		amp.x = sin(phase);
		amp.y = cos(phase);
amp.y=0;
		amp1=dot(amp,amp);        
		red=(int)(amp1*255.0);


//quantum mechanics   Schrodinger
		k=M_PI*2.0/l1;
		float position   =screen.x;
		float time 		=-screen.y/c;
		phase=position*k - time*w;
		amp.x = sin(phase);
		amp.y = cos(phase);
amp.y=0;
		amp1=dot(amp,amp);        
		green=(int)(amp1*255.0);
#endif

//composit 4d wave   (time back-fort)
		amp.x=amp.y=0;
		k=M_PI*2.0/l1a;
		phase=dot(screen,lightcone_dir1)*k;
		amp.x = sin(phase)*0.5;		amp.y = cos(phase)*0.5;

		k=M_PI*2.0/l1b;
		phase=dot(screen,lightcone_dir2)*k;//+M_PI;
		amp.x += sin(phase)*0.5;		amp.y -= cos(phase)*0.5;
//amp.y=0;
        amp1=dot(amp,amp);        
        int blue=(int)(amp1*255.0);
//green=red=0;
        pixel(x,y,(red<<16)+(green<<8)+blue);
    }

    qflush();
    getchar();
}
void elec2()
{
    float c=3e8,h=6.626e-34,m=9.1e-31,scale=1e-12,
    l1,l2,v,fi1,fi2,b,lm,f,cf,f1,f2,l1a,l2b,l1b,l2a,fia,fib,w1,w2,l4d1,l4d2;
    vec2 p1,p2,pa,pb;
    
	init_system();


	int n=100;
    scale=h/(m*c)/10.0;


    
    v=0.23*c;
    cf=c*c/v;
    fi1=atan(v/c);
    b=1.0/sqrt(1.0-v*v/(c*c));
    l1=h/(m*b*v);
    w1=m*c*c*b*2.0*M_PI/h;

   
    
    v=0.38*c;
    cf=c*c/v;
    fi2=atan(v/c);
    b=1.0/sqrt(1.0-v*v/(c*c));
    l2=h/(m*b*v);
    w2=m*c*c*b*2.0*M_PI/h;




    
    
    for(int y=0;y<550;y++)
    for(int x=0;x<800;x++)
    {
        float amp1=0,faz,l,k,w;
        vec2 amp,screen;
        
        screen.x=scale*x;
        screen.y=scale*y;


        amp.x=0;
        amp.y=0;
        for(int i=0;i<n;i++)
        {
            float t=(float)i/(n-1);

            w=w1 + (w2-w1)*t;
            l=l1 + (l2-l1)*t;

            k=M_PI*2.0/l;
            float x2   =screen.x-scale*400.0;
            float time =-screen.y/c;
            faz=x2*k - time*w;

            amp.x += sin(faz);
            amp.y += cos(faz);
        }
        amp.x/=n;//normalize
        amp.y/=n;
amp.x=0;
        amp1=dot(amp,amp);        
      int red=(int)(amp1*255.0);


        pixel(x,y,(red<<16));
    }

    qflush();
    getchar();
}
#endif

void elec3()
{
    float c=3e8,h=6.626e-34,m=9.1e-31,scale=1e-12,
	w,v,cf,fi1,gamma,l1,l4d,f,f1,f2,lm,l1a,l1b,l1c,l1d,fia,fib;
    vec2 ewave_dir,lightcone_dir1,lightcone_dir2;
    float E1,E2,f3,f4;
    
    scale=h/(m*c)/15.0;

	init_system();


    
    v=0.13*c;
    cf=c*c/v;//QM phase velocity
    fi1=atan(v/c);//angle in spacetime
    gamma=1.0/sqrt(1.0-v*v/(c*c));
    l1=h/(m*gamma*v);//De Broglie wavelength
    
    l4d=l1*sin(fi1);//spacetime  wavelength
    f=m*c*c*gamma/h;//frequency
    w=f*2.0*M_PI;//angular frequency 

#if 1
    f1=f*((c+v)/(c));//Classic Doppler
    f2=f*((c-v)/(c));//kisebb
#else
    lm=l4d/cos(fi1*2.0);    f=sqrt(c*c+v*v)/lm;//        f=f*(cf/(cf-v));
    f1=f*(c/(c-v));//Classic Doppler
    f2=f*(c/(c+v));//kisebb
#endif
//    f=(f1+f2)/2.0;

    l1a=c/f1*sin(45.0*iradian);//a component on lightcone
    l1b=c/f2*sin(45.0*iradian);//b component on lightcone
    
 


    
    v=0.61*c;
    cf=c*c/v;//QM phase velocity
    fi1=atan(v/c);//angle in spacetime
    gamma=1.0/sqrt(1.0-v*v/(c*c));
    l1=h/(m*gamma*v);//De Broglie wavelength
    
    l4d=l1*sin(fi1);//spacetime  wavelength
    f=m*c*c*gamma/h;//frequency
    w=f*2.0*M_PI;//angular frequency 

#if 1
    f1=f*((c+v)/(c));//Classic Doppler
    f2=f*((c-v)/(c));//kisebb
#else
    lm=l4d/cos(fi1*2.0);    f=sqrt(c*c+v*v)/lm;//        f=f*(cf/(cf-v));
    f1=f*(c/(c-v));//Classic Doppler
    f2=f*(c/(c+v));//kisebb
#endif
//    f=(f1+f2)/2.0;

    l1c=c/f1*sin(45.0*iradian);//a component on lightcone
    l1d=c/f2*sin(45.0*iradian);//b component on lightcone
    

    
    fia=45.0*iradian;
    lightcone_dir1.x=sin(fia);
    lightcone_dir1.y=cos(fia);
    fib=-45.0*iradian;
    lightcone_dir2.x=sin(fib);
    lightcone_dir2.y=cos(fib);

    
    
    for(int y=0;y<550;y++)
    for(int x=0;x<800;x++)
    {
	    int red=0,green=0;
		float amp1=0,phase,l,k;
		vec2 amp,screen;
        
		screen.x=scale*x;
		screen.y=scale*y;

w=0.25;//0.25
//composit 4d wave   (time back-fort)
		amp.x=amp.y=0;
		int n=20;
     //   for(int i=0;i<n;i++)
        {
//            float t=(float)i/(n-1);
//            t=1.0+0.1*t;

//	float t=1.0;	if(y<250) t=0.0;

#if 1		
		k=M_PI*2.0/l1a;
		phase=dot(screen,lightcone_dir1)*k;
		amp.x += sin(phase)*w;		amp.y += cos(phase)*w;

		k=M_PI*2.0/l1b;
		phase=dot(screen,lightcone_dir2)*k;//+M_PI;
		amp.x += sin(phase)*w;		amp.y -= cos(phase)*w;
#endif
#if 0
		vec2 amp2;
		k=M_PI*2.0/l1c;
		phase=dot(screen,lightcone_dir1)*k;
	//	amp.x += sin(phase)*w;		amp.y += cos(phase)*w;

		k=M_PI*2.0/l1d;
		phase=dot(screen,lightcone_dir2)*k;//+M_PI;
		amp.x += sin(phase)*w;		amp.y -= cos(phase)*w;
#endif		
}
		//amp.x/=n;		amp.y/=n;
		
//amp.y=0;
        amp1=dot(amp,amp);        
        int blue=(int)(amp1*255.0);

        pixel(x,y,(red<<16)+(green<<8)+blue);
    }

    for(int y=0;y<550;y++) pixel(y,y,0x00ff00);

    qflush();
    getchar();
}
//__________________________________________________________________


//__________________________________________________________________

double c=3e8,x3,t3=-1.0,x4,t4=-1.0;
double wide=30.0;



void setorigin(int &x,int &y)
{
	x+=300;
	y=630-y;    
}
void line2(int x1,int y1,int x2,int y2,int color)
{
	setorigin(x1,y1);
	setorigin(x2,y2);
	line(x1,y1,x2,y2,color);
}


void lorentztransformation(double x1,double t1,double &x2,double &t2,double v)
{
	double gamma=1.0/sqrt(1.0-v*v/(c*c));

	x2=(x1-v*t1)*gamma;
	t2=(t1-v*x1/(c*c))*gamma;
}	
void drawworldline(double x1,double t1,double v)
{
	double x2,t2;
	lorentztransformation(x1,t1,x2,t2,v);

	if(t4!=-1.0)	line2(floor(x2),floor(t2*c),floor(x4),floor(t4*c),0x880000);
	if(t3!=-1.0)	line2(floor(x2),floor(t2*c),floor(x3),floor(t3*c),0xffff00);
	//BUG would be nice if yellow line be on top always

	x4=x3;t4=t3;
	x3=x2;t3=t2;//store last
}
void drawtimecoords(double v,int color)
{
	double x2,t2,x2b,t2b;
	double x=0,time=0;
	
	for(int i=0;i<16;i++)
	{
		lorentztransformation(x-100,time,x2,t2,v);
		lorentztransformation(x+400,time,x2b,t2b,v);
		
		line2(x2,t2*c,x2b,t2b*c,color);
		time+=wide/c;
	}
}
void drawlightclock(double xposition,double time,double v,int n=16)
{
	t3=-1.0;t4=-1.0;//reset
	int i=0;
	
    while(i<=n)
    {
    	drawworldline(xposition,time,v);
    	xposition+=wide; //v*dt  v=c  c*wide/c=wide
    	time+=wide/c;//dt=wide/c
    	i++;

		if(i<=n)
		{
	    	drawworldline(xposition,time,v);
	    	xposition-=wide; 
	    	time+=wide/c;
	    }
	    i++;
    }
}
int lclock()
{
	init_system();

	double v=-0.6*c;//60% of lightspeed
//v=0.6*c;	
	drawtimecoords(0.0,0x005500);
	drawtimecoords(v,0x000055);
	
	drawlightclock(0.0,0.0,0.0);
	drawlightclock(200.0,0.0,0.0);
	drawlightclock(0.0,0.0,v);	
	drawlightclock(200.0,0.0,v);//0.0000004


    qflush();
    getchar();
        
    return 0;
}


int lclock2()
{
	init_system();

	double v=-0.6*c;//60% of lightspeed
//v=0.6*c;	
	drawtimecoords(0.0,0x005500);
//	drawtimecoords(v,0x000055);
	
	drawlightclock(0.0,0.0,0.0,20);
	drawlightclock(0.0,0.0,v,8);
	double x5,t5;	
	lorentztransformation(x3,t3,x5,t5,v);
	drawlightclock(x5,t5,-v,8);	

double t=1.0/sqrt(1.0-v*v/(c*c));
printf("moving clock in external view %f \n",t*8*2);

    qflush();
    getchar();
        
    return 0;
}

int lclock3()
{
	init_system();

	wide=25.0;//30 too large
	
	double v=-0.6*c;//60% of lightspeed
//v=0.6*c;	
	drawtimecoords(0.0,0x005500);
//	drawtimecoords(v,0x000055);
	
	drawlightclock(0.0,0.0,v,20);

	drawlightclock(0.0,0.0,0.0,8);
	double u=v;
	double v2=(v+u)/(1.0+v*u/(c*c));//v+v!   https://en.wikipedia.org/wiki/Velocity-addition_formula
	double x5,t5;		
	lorentztransformation(x3,t3,x5,t5,-v2);
	drawlightclock(x5,t5,v2,8);	



double t=1.0/sqrt(1.0-v*v/(c*c));
printf("moving clock in external view %f \n",t*8*2);

    qflush();
    getchar();
        
    return 0;
}


//__________________________________________________________________

//__________________________________________________________________


int main()
{
/*
    float c=3e8,h=6.626e-34,m=9.1e-31,scale=1e-12,
	w,v,cf,fi1,gamma,l1,l4d,f,f1,f2,lm,l1a,l1b,fia,fib;
    vec2 ewave_dir,lightcone_dir1,lightcone_dir2;

    v=0.58*c;
    v=0.87*c;
    cf=c*c/v;//QM phase velocity
    fi1=atan(v/c);//angle in spacetime
    gamma=1.0/sqrt(1.0-v*v/(c*c));
    l1=h/(m*gamma*v);//De Broglie wavelength
    
    l4d=l1*sin(fi1);//spacetime  wavelength
printf("%e \n",l4d);

l1=h/(m*c)*cos(67.5*iradian);//l1=h/(m*c)*c/sqrt(c*c+c*c);
printf("%e \n",l1);
*/

dslitepr();
//dslitepr2();
//qwp();
//epr();
//elec();
//elec3();
//lclock();
//lclock2();
//lclock3();

  qflush();
    getchar();
  

}

/* 

        w=2*m*g/(r*c*c*4)/2;
        c2=c*(1-w)/((1+w)*(1+w)*(1+w));


 
*/
 
#if 0


#include <stdio.h>
#include <math.h>
#include <X11/Xlib.h>


Display *dpy;
Window win;
GC gc;


typedef long double float1;



struct float3
{
    float1 x,y,z;
    
    float3() {x=0;y=0;z=0;};
    float3(float1 x1,float1 y1,float1 z1) {x=x1;y=y1;z=z1;};
    float3 operator + (float3 v) {float3 v2;v2.x = x+v.x;v2.y = y+v.y;v2.z = z+v.z;return v2;};
    float3 operator - (float3 v) {float3 v2;v2.x = x-v.x;v2.y = y-v.y;v2.z = z-v.z;return v2;};
    float3 operator * (float1 s) {float3 v2;v2.x = x*s;v2.y = y*s;v2.z = z*s;return v2;};
};
float1 dot(float3 v1,float3 v2) {    float1 t = v1.x*v2.x+v1.y*v2.y+v1.z*v2.z;return t;}
float1 length(float3 v1) {    return sqrtl(dot(v1,v1));}
float3 normalize(float3 v1){    return v1*(1.0/length(v1));}





float1 rs,dt,rs2,g,c,m,r1,v1,
    c2 = 1.0,
    skala2 = 1.0,
    skala;
float3 pozicio, sebesseg;
int n;



void pont(float3 v1,int szin)
{
    v1.x*=skala2;
    v1.y*=skala2;
    v1.z*=skala2;
    
    XSetForeground(dpy,gc,szin);
    XDrawPoint(dpy, win, gc, 500+v1.x,300+v1.y);
}

void schwarzschild()
{
    float1 r,p,t,T, dr,dp, ddr,ddT,ddp,  dT;
    
    skala = c;//  c=1
    c2  = c /skala;
    rs2 = rs/skala;

    pozicio.x  = pozicio.x/skala;
    pozicio.y  = pozicio.y/skala;
    pozicio.z  = pozicio.z/skala;
    sebesseg.x  = sebesseg.x/skala;
    sebesseg.y  = sebesseg.y/skala;
    sebesseg.z  = sebesseg.z/skala;



    r  = sqrtl(pozicio.x*pozicio.x + pozicio.y*pozicio.y);
    p  = atan2(pozicio.y,pozicio.x);
    t = 0.0;

    dr = (pozicio.x * sebesseg.x + pozicio.y * sebesseg.y) / r;
    dp = (pozicio.x * sebesseg.y - pozicio.y * sebesseg.x) / (r*r);
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
        

        pozicio.x = r*cos(p)*skala;
        pozicio.y = r*sin(p)*skala;

        
        pont(pozicio,0x00ff00);
    }
}

void newton()
{
    for(int j = 0;j<n;j++)
    {
        pozicio.x = pozicio.x + sebesseg.x*dt;
        pozicio.y = pozicio.y + sebesseg.y*dt;
        pozicio.z = pozicio.z + sebesseg.z*dt;

        float3 av = pozicio;
        float1 r=length(av);
        av.x/=r;
        av.y/=r;
        av.z/=r;
        
        float1 a = -m*g/(r*r);
        sebesseg.x = sebesseg.x + av.x*a*dt;
        sebesseg.y = sebesseg.y + av.y*a*dt;
        sebesseg.z = sebesseg.z + av.z*a*dt;

        pont(pozicio,0x0000ff);
    }
}

void refract()
{
    float1 r,c1,c2, cx,cy,w;
    float3 pozicio4,sebesseg4;//4d  x,y,time   z nincs
    

      pozicio.z=0;
      sebesseg.z=0;

    r=length(pozicio)-rs;
    w=2*m*g/(r*c*c*4)/2;
    c1=c*(1-w)/((1+w)*(1+w)*(1+w));

 
    pozicio4=pozicio;
    sebesseg4=sebesseg;
    sebesseg4.z=sqrtl(c1*c1 - sebesseg.x*sebesseg.x - sebesseg.y*sebesseg.y);// .z == time
    
    for(int j = 0;j<n;j++)
    {
        pozicio4 = pozicio4+sebesseg4*dt;
        pozicio=pozicio4;
        pozicio.z=0;

        float3 N = normalize(pozicio);
        float3 T = normalize(sebesseg4);          
        T=normalize(T-N*dot(N,T));//  T meroleges N re

        r=length(pozicio)-rs;
        w=2*m*g/(r*c*c*4)/2;
        c2=c*(1-w)/((1+w)*(1+w)*(1+w));

      
        cx=dot(T,sebesseg4);
        cy=dot(N,sebesseg4);
        float1 cy2=cy;

        cx=cx*c2*c2/(c1*c1);
        cy=(c2*c2-cx*cx);
        if(cy<0.0) cy=-sqrtl(-cy);
        else       cy= sqrtl(cy);
        if(cy2<0.0) cy=-cy;

        sebesseg4=normalize(T*cx + N*cy)*c2;
        c1=c2;


        pont(pozicio,0xff0000);
    }
}

void setup()
{
     g = 6.674e-11;
     c = 2.997e8;
     c2 =1.0;
     m = 1.98e30;

    rs = 2.0*m*g/(c*c);
    r1 = rs*500.0;    dt = 2e-5;    n = 2000000;
//    r1 = rs*50.0;    dt = 2e-7;    n = 2000000;
//    r1 = rs*15.0;    dt = 2e-8;    n = 4000000;
    v1 = sqrt(m*g/r1)*1.2;

    skala2 = 100.0/r1;
    pozicio = float3(r1,0.0f,0.0f);
   
    float1 alfa=61.0*M_PI/180.0;
    sebesseg = float3(v1*cos(alfa),v1*sin(alfa),0.0f);
}
int main()
{
    dpy  =  XOpenDisplay((0));
    win  =  XCreateSimpleWindow(dpy, DefaultRootWindow(dpy), 0,0, 1200, 800, 0,0,0);
   
    XSelectInput(dpy, win, StructureNotifyMask);
    XMapWindow(dpy, win);
   
    gc  =  XCreateGC(dpy, win, 0, (0));
    XSetForeground(dpy,gc,0);
   
    for(;;) {    XEvent e;    XNextEvent(dpy, &e);    if (e.type  ==  MapNotify)break;    }
   

    
    setup();newton();   
    setup();schwarzschild();
    setup();refract();
    
    
    XFlush(dpy);
    getchar();

    return 0;

}





#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <X11/Xlib.h>



Display *dpy;
Window win;
GC gc;


void pixel(int x,int y,int color)
{
    XSetForeground(dpy,gc,color);
    XDrawPoint(dpy, win, gc, x,y);
}
void line(int x1,int y1,int x2,int y2,int color)
{
    XSetForeground(dpy,gc,color);
    XDrawLine(dpy, win, gc, x1,y1,x2,y2);
}
void init_system()
{
    dpy = XOpenDisplay(0);
    win = XCreateSimpleWindow(dpy, DefaultRootWindow(dpy), 0,0, 900, 650, 0,0,0);
    
    XSelectInput(dpy, win, StructureNotifyMask);
    XMapWindow(dpy, win);
    gc = XCreateGC(dpy, win, 0, 0);
    
    for(;;) { XEvent e; XNextEvent(dpy, &e); if (e.type == MapNotify) break; }
}





double cos2(double a) 
{
	a=cos(a);
	return a*a;
}
double sin2(double a) 
{
	a=sin(a);
	return a*a;
}


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



int bell()
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
    	S=S2qm(a,b,a2,b2 );    			printf("qm: %f \n",S); 			pixel(angle*4+1,200-S*50.0,0xff0000);//red   
	    XFlush(dpy);
	}

    XFlush(dpy);
    getchar();
}

int main()
{
	bell();
}




//g++ phys2.cpp -O3 -lX11

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <X11/Xlib.h>

#include <gmp.h>
#include <gmpxx.h>

typedef mpf_class f32;




Display *dpy;
Window win;
GC gc;


void pixel(int x,int y,int color)
{
    XSetForeground(dpy,gc,color);
    XDrawPoint(dpy, win, gc, x,y);
}
void line(int x1,int y1,int x2,int y2,int color)
{
    XSetForeground(dpy,gc,color);
    XDrawLine(dpy, win, gc, x1,y1,x2,y2);
}
void init_system()
{
    dpy = XOpenDisplay(0);
    win = XCreateSimpleWindow(dpy, DefaultRootWindow(dpy), 0,0, 900, 650, 0,0,0);
    
    XSelectInput(dpy, win, StructureNotifyMask);
    XMapWindow(dpy, win);
    gc = XCreateGC(dpy, win, 0, 0);
    
    for(;;) { XEvent e; XNextEvent(dpy, &e); if (e.type == MapNotify) break; }
}





double cos2(double a) 
{
	a=cos(a);
	return a*a;
}
double cos3(double a) 
{
	a=cos(a);
	return fabs(a);
///	return a;
}
double sin2(double a) 
{
	a=sin(a);
	return a*a;
}
double sin3(double a) 
{
	a=sin(a);
	return fabs(a);
//	return a;
}







f32 h,R,c=2.9979e8,ly,G=6.7e-11,m1=2e30,m2=6e24,r,q,one=1.0,two=2.0,
		f,l,z,mpc,w,v,e,k,x,F,D,m,E,T,K,rs1,rs2,v1,v2,A,n,V,dr,dv1,dv2,dp,dt,r2,r3,dt2,h2;


f32 sqr(f32 n)
{
	return n*n;
}
f32 force(int n)
{
	if(n==0) return (e*e*k);
	if(n==1) return (e*e*k*137.0);
	if(n==2) return (m*m*G);
	return f32(0);
}

int main()
{
	mpf_set_default_prec(512*32);
//	bell();

	r=150e9;
	ly=c*24.0*3600.0*365.0;
	mpc=ly*3.2e6;

	R=1.0/4.0*ly;//R*=1e5;
m2=m1*29.0;
m1=m1*36.0;
r=210e3;//r/=6.28;	printf("r  %.12e \n",r.get_d());
R=mpc*410.0;


	h=-1.0/R * G*G/(c*c*c*c) * 4.0*m1*m2/r;
	h= G*G/(c*c*c*c) * 4.0*m1*m2/r;
	printf("h  %.12e \n",h.get_d());

	rs1= (m1*2.0*G/(c*c));
	rs2= (m2*2.0*G/(c*c));//	printf("rs %.12e \n",rs1.get_d());
	v1=sqrt(G*m1/r);
	v2=sqrt(G*m2/r);



	
	
//suruseg
#if 0
q=f32(2e9)/(c*c);   //v=sqrt(2e9/ro)  r=2e9/v2
	printf("%.15e \n\n",q.get_d());

//v	
	q=1000.0;//h2o
	v=sqrt(f32(2e9)/q);
	printf("v %.15e \n",v.get_d());
	q=13500.0;//hg
	v=sqrt(f32(30e9)/q);
	printf("v %.15e \n",v.get_d());
	q=7860.0;//Fe
	v=sqrt(f32(205e9)/q);//K
	printf("v? %.15e \n",v.get_d());//5100

//ampl
	f=880.0;
	w=2.0*M_PI*f;
	v=340.0;
	q=v/w;
	printf("%.15e \n\n",q.get_d());
#endif

/*
	f=220.0;
	w=2.0*M_PI*f;
	v=c;
	q=v/w;
	printf("%.15e \n",q.get_d());
	q/=R;//5*
	printf("%.15e \n\n",q.get_d());
	*/


	int u=0;
//spring
	m=9.1e-31; //15.8 GeV  +strong     2147GeV elec
//m*=105.0*2.0;//3333 GeV +strong     456.7 GeV
//m=1e-7;//GR  8e22GeV
//m=1e-41;

m/=1e15;//3e4;
//m1*=1e4;


	e=1.6e-19;
	k=9e9;
	v=c;
			//	x=52.9e-12/137.0/137.0*4.0;//	printf("%.15e \n",x.get_d());
//	x=4.0e-14; //STRONG
//x=0.3e-15; //Pr  2000 GeV
//x=1e-35;//GR



	//F=e*e*k/(x*x)*137.0;//STRONG 18   15GeV
//F=e*e*k/(x*x);//Pr  23 
//F=m*m*G/(x*x);//GR  105   5e22GeV lol
//	D=F/x;	f=f32(1.0/(2.0*M_PI))*sqrt(D/m);	T=one/f;	v=x/(T/2.0);
//v=x*2.0*sqrt(e*e*k*137.0/(x*x*x*m))/(2.0*M_PI) ;      //v=2xsqrt(eek137/(xxxm))/(2pi)
	x=force(u)/(sqr(c*M_PI)*m);  //!!!!!!!111
	
	printf("x %.15e \n",x.get_d());
v=x*sqrt(force(u)/(x*x*x*m))/(M_PI) ;	
	printf("v %.15e \n",v.get_d());
	q=one/x;	q=q*q*q; q*=m; K=v*v*q/1e9;//K
	printf("k %.15e \n\n",K.get_d());
	

	
	F=force(u)/(x*x);	D=F/x; f=f32(1.0/(2.0*M_PI))*sqrt(D/m);//re
//	printf("%.15e \n",f.get_d());
	h=6.626e-34;
	E=f*h/(e*1e6);//MeV
	printf("E %.15e \n",E.get_d());

	E=D*x*x*0.5;
//	printf("E %.15e \n",E.get_d());




	
#if 0
//FE!
	m=9.1e-31*1836.0*55.8;//Fe
	e=1.6e-19;
	k=9e9;
	x=228.05e-12;//	x=52.9e-12*4.6;
	printf("x %.15e \n",x.get_d());
	F=e*e*k/(x*x)*26.0*0.917;//92
	D=F/x;	f=f32(1.0/(2.0*M_PI))*sqrt(D/m);	T=one/f;	v=x/(T/2.0);
	printf("v %.15e \n",v.get_d());
	q=one/x;	q=q*q*q; q*=m; //q=7860
	K=v*v*q/1e9;//K  v=sqrt(K/q)
	printf("K %.15e \n",K.get_d());
q=7860.0;
q=q/m;
q=one/q;
	printf("x %.15e \n",pow(q.get_d(),1.0/3.0));//x
#endif

	
}




//g++ phys2.cpp -O3 -lX11

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <X11/Xlib.h>

#include <gmp.h>
#include <gmpxx.h>

typedef mpf_class f32;




Display *dpy;
Window win;
GC gc;


void pixel(int x,int y,int color)
{
    XSetForeground(dpy,gc,color);
    XDrawPoint(dpy, win, gc, x,y);
}
void line(int x1,int y1,int x2,int y2,int color)
{
    XSetForeground(dpy,gc,color);
    XDrawLine(dpy, win, gc, x1,y1,x2,y2);
}
void init_system()
{
    dpy = XOpenDisplay(0);
    win = XCreateSimpleWindow(dpy, DefaultRootWindow(dpy), 0,0, 900, 650, 0,0,0);
    
    XSelectInput(dpy, win, StructureNotifyMask);
    XMapWindow(dpy, win);
    gc = XCreateGC(dpy, win, 0, 0);
    
    for(;;) { XEvent e; XNextEvent(dpy, &e); if (e.type == MapNotify) break; }
}





double cos2(double a) 
{
	a=cos(a);
	return a*a;
}
double cos3(double a) 
{
	a=cos(a);
	return fabs(a);
///	return a;
}
double sin2(double a) 
{
	a=sin(a);
	return a*a;
}
double sin3(double a) 
{
	a=sin(a);
	return fabs(a);
//	return a;
}





f32 h,R,c=2.9979e8,ly,G=6.7e-11,m1=2e30,m2=6e24,r,q,one=1.0,two=2.0,
		f,l,z,mpc,w,v,e,k,x,F,D,m,E,T,K,rs1,rs2,v1,v2,A,n,V,dr,dv1,dv2,dp,dt,r2,r3,dt2,h2;


f32 sqr(f32 n)
{
	return n*n;
}
f32 force(int n)
{
	if(n==0) return (e*e*k);
	if(n==1) return (e*e*k*137.0);
	if(n==2) return (m*m*G);
	return f32(0);
}

int main()
{
	mpf_set_default_prec(512*32);
//	bell();

	r=150e9;
	ly=c*24.0*3600.0*365.0;
	mpc=ly*3.2e6;

	R=1.0/4.0*ly;//R*=1e5;
m2=m1*29.0;
m1=m1*36.0;
r=210e3;//r/=6.28;	printf("r  %.12e \n",r.get_d());
R=mpc*410.0;


	h=-1.0/R * G*G/(c*c*c*c) * 4.0*m1*m2/r;
	h= G*G/(c*c*c*c) * 4.0*m1*m2/r;
	printf("h  %.12e \n",h.get_d());

	rs1= (m1*2.0*G/(c*c));
	rs2= (m2*2.0*G/(c*c));//	printf("rs %.12e \n",rs1.get_d());
	v1=sqrt(G*m1/r);
	v2=sqrt(G*m2/r);



	
	
//suruseg
#if 0
q=f32(2e9)/(c*c);   //v=sqrt(2e9/ro)  r=2e9/v2
	printf("%.15e \n\n",q.get_d());

//v	
	q=1000.0;//h2o
	v=sqrt(f32(2e9)/q);
	printf("v %.15e \n",v.get_d());
	q=13500.0;//hg
	v=sqrt(f32(30e9)/q);
	printf("v %.15e \n",v.get_d());
	q=7860.0;//Fe
	v=sqrt(f32(205e9)/q);//K
	printf("v? %.15e \n",v.get_d());//5100

//ampl
	f=880.0;
	w=2.0*M_PI*f;
	v=340.0;
	q=v/w;
	printf("%.15e \n\n",q.get_d());
#endif

/*
	f=220.0;
	w=2.0*M_PI*f;
	v=c;
	q=v/w;
	printf("%.15e \n",q.get_d());
	q/=R;//5*
	printf("%.15e \n\n",q.get_d());
	*/


	int u=0;
//spring
	m=9.1e-31; //15.8 GeV  +strong     2147GeV elec
//m*=105.0*2.0;//3333 GeV +strong     456.7 GeV
//m=1e-7;//GR  8e22GeV
//m=1e-41;

m/=1e15;//3e4;
//m1*=1e4;


	e=1.6e-19;
	k=9e9;
	v=c;
			//	x=52.9e-12/137.0/137.0*4.0;//	printf("%.15e \n",x.get_d());
//	x=4.0e-14; //STRONG
//x=0.3e-15; //Pr  2000 GeV
//x=1e-35;//GR



	//F=e*e*k/(x*x)*137.0;//STRONG 18   15GeV
//F=e*e*k/(x*x);//Pr  23 
//F=m*m*G/(x*x);//GR  105   5e22GeV lol
//	D=F/x;	f=f32(1.0/(2.0*M_PI))*sqrt(D/m);	T=one/f;	v=x/(T/2.0);
//v=x*2.0*sqrt(e*e*k*137.0/(x*x*x*m))/(2.0*M_PI) ;      //v=2xsqrt(eek137/(xxxm))/(2pi)
	x=force(u)/(sqr(c*M_PI)*m);  //!!!!!!!111
	
	printf("x %.15e \n",x.get_d());
v=x*sqrt(force(u)/(x*x*x*m))/(M_PI) ;	
	printf("v %.15e \n",v.get_d());
	q=one/x;	q=q*q*q; q*=m; K=v*v*q/1e9;//K
	printf("k %.15e \n\n",K.get_d());
	

	
	F=force(u)/(x*x);	D=F/x; f=f32(1.0/(2.0*M_PI))*sqrt(D/m);//re
//	printf("%.15e \n",f.get_d());
	h=6.626e-34;
	E=f*h/(e*1e6);//MeV
	printf("E %.15e \n",E.get_d());

	E=D*x*x*0.5;
//	printf("E %.15e \n",E.get_d());




	
#if 0
//FE!
	m=9.1e-31*1836.0*55.8;//Fe
	e=1.6e-19;
	k=9e9;
	x=228.05e-12;//	x=52.9e-12*4.6;
	printf("x %.15e \n",x.get_d());
	F=e*e*k/(x*x)*26.0*0.917;//92
	D=F/x;	f=f32(1.0/(2.0*M_PI))*sqrt(D/m);	T=one/f;	v=x/(T/2.0);
	printf("v %.15e \n",v.get_d());
	q=one/x;	q=q*q*q; q*=m; //q=7860
	K=v*v*q/1e9;//K  v=sqrt(K/q)
	printf("K %.15e \n",K.get_d());
q=7860.0;
q=q/m;
q=one/q;
	printf("x %.15e \n",pow(q.get_d(),1.0/3.0));//x
#endif

	
}






int lorentzcontraction()
{
	init_system();



	for(int i=0;i<500;i++)
	{
	v0=c*i/500.0;
	u=c/137.0;
	
	v=u;
	gamma=1/sqrt(1.0-v*v/(c*c));
	L0=h/(m*v*gamma*2.0*M_PI);
	//printf("%e \n",L0);



	v=v0;
	v=(v+u)/(1.0+v*u/(c*c));
	gamma=1/sqrt(1.0-v*v/(c*c));
	L1=h/(m*v*gamma*2.0*M_PI);

	v=v0;
	u=-c/137.0;
	v=(v+u)/(1.0+v*u/(c*c));
	gamma=1/sqrt(1.0-v*v/(c*c));
	L2=h/(m*v*gamma*2.0*M_PI);
	f=(c/L1-c/L2)/2.0;
//	L=2.0*(L2*L1)/(L2-L1);//wavepacket of electron of atom
	L=c/f;
	L5=L;
//	printf("%.12e ",L);
	pixel(i,400-L*5e12,0x00ff00);


	v=v0;
	gamma=1/sqrt(1.0-v*v/(c*c));
	L=L0/gamma;//Lorentz contraction
	printf("%.12e \n",L/L5);
	
	//shifted in y direction by 1 pixel
	pixel(i,401-L*5e12,0xff0000);
	}
  //  XFlush(dpy);    getchar();
        
    return 0;
}
int QMwavedecomposition()
{
	init_system();

	double c=3e8,v,gamma,h=6.626e-34,L,m=9.1e-31,f0,f1,f2,f,L0,L2;
	L0=h/(m*c*M_PI);//Compton wavelength
	
	
	for(int i=0;i<500;i++)
	{
		v=c*i/500.0;
		gamma=1/sqrt(1.0-v*v/(c*c));
		L2=h/(m*v*gamma*2.0*M_PI); //QM   De Broglie wavelength
		pixel(i,400-L2*25e12,0xff0000);

//QM wave decomposition
		f0=c/L0/gamma;
	    f1=f0*(c/(c-v));//Classic Doppler   /wave source is moving/
    	f2=f0*(c/(c+v));
	    f=(f1-f2);//wavepacket  2 waves in diffecent directions
		L=c/f;
		pixel(i,401-L*25e12,0x00ff00);
//printf("%e \n",L/L2);		
	}

    XFlush(dpy);
    getchar();
        
    return 0;
}
int main()
{
//	QMwavedecomposition();
	lorentzcontraction();

}


#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 


struct complex 
{ 
double real , img; 
complex() {real = img = 0;}; 
complex(double r, double i) {real =r ; img = i;}; 
void set(double r, double i) {real=r ; img = i;}; 
void conjugate() {img *= -1;} 
complex operator + (complex c) 
{ 
	complex e;
	e.real = real + c.real; 
	e.img  = img  + c.img; 
	return e; 
} 
complex operator - (complex c) 
{ 
	complex e;
	e.real = real - c.real; 
	e.img  = img  - c.img ; 

	return e; 
} 
complex operator * (complex c) 
{ 
	complex e;
	e.real = real*c.real - img *c.img; 
	e.img  = img *c.real + real*c.img; 

	return e; 
} 
complex operator / (complex c) 
{ 
	complex e;
	double denominator = c.real*c.real + c.img*c.img; 

	e.real = (real * c.real + img * c.img) / denominator; 
	e.img = (img * c.real - real *c.img ) / denominator; 
	return e; 
} 
}; 
double alpha, p, m_x, m_y, m_z, n_x, n_y, n_z; 
complex sigma_z[2][2], sigma_x[2][2], sigma_y[2][2]; 
complex sigmaDOTm[2][2], sigmaDOTn[2][2], spin_n[2], spin_m[2]; 
complex v1[2], v2[2], v3[2], n[3], m[3], P, P2, normalization_factor; 


double probability (complex * c1, complex *c2 ) 
{ 
	complex c3[2], P; 

	c3[0] = c1[0]; 
	c3[1] = c1[1]; 
	c3[0].conjugate(); // c3 =  <c1| 
	c3[1].conjugate(); 

	P = c3[0]*c2[0] + c3[1]*c2[1]; // <c1|c2> 
	P = P*P; // ^2  

	return P.real; 
} 

void transf()
{
//sigma dot n and m
	sigmaDOTn[0][0] = sigma_x[0][0]*n[0] + sigma_y[0][0]*n[1] + sigma_z[0][0]*n[2];   
	sigmaDOTn[0][1] = sigma_x[0][1]*n[0] + sigma_y[0][1]*n[1] + sigma_z[0][1]*n[2];   
	sigmaDOTn[1][0] = sigma_x[1][0]*n[0] + sigma_y[1][0]*n[1] + sigma_z[1][0]*n[2];   
	sigmaDOTn[1][1] = sigma_x[1][1]*n[0] + sigma_y[1][1]*n[1] + sigma_z[1][1]*n[2];   
	printf( "? >%.2f/%.2f  ", sigmaDOTn[0][0].real, sigmaDOTn[0][0].img); 
	printf( "? >%.2f/%.2f  \n", sigmaDOTn[0][1].real, sigmaDOTn[0][1].img); 
	printf( "? >%.2f/%.2f  ", sigmaDOTn[1][0].real, sigmaDOTn[1][0].img); 
	printf( "? >%.2f/%.2f  \n", sigmaDOTn[1][1].real, sigmaDOTn[1][1].img); 

	v2[0] = sigmaDOTn[0][0]*spin_n[0] + sigmaDOTn[0][1]*spin_n[1]; 
	v2[1] = sigmaDOTn[1][0]*spin_n[0] + sigmaDOTn[1][1]*spin_n[1]; 
	printf( "? sigmaDOTn>%.2f/%.2f   %.2f/%.2f \n", spin_n[0].real, spin_n[0].img, spin_n[1].real, spin_n[1].img); 
	printf( "? sigmaDOTn>%.2f/%.2f   %.2f/%.2f \n", v2[0].real, v2[0].img, v2[1].real, v2[1].img); 

	spin_n[0]=v2[0];
	spin_n[1]=v2[1];
}



int main () 
{ 
//Pauli matrices
	sigma_x[0][0].set(0,0); sigma_x[0][1].set(1,0); 
	sigma_x[1][0].set(1,0); sigma_x[1][1].set(0,0); 

	sigma_y[0][0].set(0,0); sigma_y[0][1].set(0, -1); 
	sigma_y[1][0].set(0,1); sigma_y[1][1].set(0,0); 

	sigma_z[0][0].set(1,0); sigma_z[0][1].set(0,0); 
	sigma_z[1][0].set(0,0); sigma_z[1][1].set(-1,0); 


//3d vectors
	alpha = 77.0* M_PI/180.0; 
	m_x = 0; 
	m_y = sin(alpha); 
	m_z = cos(alpha); 

	m[0].set(m_x, 0); 
	m[1].set(m_y, 0); 
	m[2].set(m_z, 0); 
	
	n_x = 0; 
	n_y = 0; 
	n_z = 1; 
	n[0].set(n_x, 0); 
	n[1].set(n_y, 0); 
	n[2].set(n_z, 0); 

	sigmaDOTn[0][0] = sigma_x[0][0]*n[0] + sigma_y[0][0]*n[1] + sigma_z[0][0]*n[2];   
	sigmaDOTn[0][1] = sigma_x[0][1]*n[0] + sigma_y[0][1]*n[1] + sigma_z[0][1]*n[2];   
	sigmaDOTn[1][0] = sigma_x[1][0]*n[0] + sigma_y[1][0]*n[1] + sigma_z[1][0]*n[2];   
	sigmaDOTn[1][1] = sigma_x[1][1]*n[0] + sigma_y[1][1]*n[1] + sigma_z[1][1]*n[2];   
	sigmaDOTm[0][0] = sigma_x[0][0]*m[0] + sigma_y[0][0]*m[1] + sigma_z[0][0]*m[2];  
	sigmaDOTm[0][1] = sigma_x[0][1]*m[0] + sigma_y[0][1]*m[1] + sigma_z[0][1]*m[2];  
	sigmaDOTm[1][0] = sigma_x[1][0]*m[0] + sigma_y[1][0]*m[1] + sigma_z[1][0]*m[2];  
	sigmaDOTm[1][1] = sigma_x[1][1]*m[0] + sigma_y[1][1]*m[1] + sigma_z[1][1]*m[2];  




//spins
	spin_n[0].set(1,0); 
	spin_n[1].set(0,0); 

	n[0].set(0, 0); 
	n[1].set(0, 0); 
	n[2].set(1, 0); 
	transf();

	n[0].set(0, 0); 
	n[1].set(0, 0); 
	n[2].set(-1, 0); 
	transf();

	n[0].set(0, 0); 
	n[1].set(0, 0); 
	n[2].set(1, 0); 
	transf();

	n[0].set(0, 0); 
	n[1].set(0, 0); 
	n[2].set(-1, 0); 
	transf();

	n[0].set(0, 0); 
	n[1].set(0, 0); 
	n[2].set(1, 0); 
	transf();




	v2[0] = sigmaDOTn[0][0]*spin_n[0] + sigmaDOTn[0][1]*spin_n[1]; 
	v2[1] = sigmaDOTn[1][0]*spin_n[0] + sigmaDOTn[1][1]*spin_n[1]; 
	printf( "| sigmaDOTn>%e %e %e %e \n", v2[0].real, v2[0].img, v2[1].real, v2[1].img); 

	

#if 0
p = (1.0 + cos(alpha)) / 2.0; //check reverse
p = sqrt(p); 
spin_m[0].set(p, 0); 
spin_m[1].set(sqrt (1 - p*p), 0); 
printf( "spin         %e %e  \n", spin_m[0].real, spin_m[1].real); 
#endif 

	spin_m[0].set(1,0); 
	spin_m[1] = complex(1 - m_z, 0) / complex(m_x, -m_y); // P = ((1 - m_z) / (m_x, m_y * i)); 
	normalization_factor.set(sqrt(((1 + m_z) / 2)), 0); 
	spin_m[0] = spin_m[0]*normalization_factor; 
	spin_m[1] = spin_m[1]*normalization_factor; 
	printf( "spin         %e %e   \n", spin_m[0].real, spin_m[1].real); 



	v1[0] = sigmaDOTm[0][0]*spin_m[0] + sigmaDOTm[0][1]*spin_m[1]; 
	v1[1] = sigmaDOTm[1][0]*spin_m[0] + sigmaDOTm[1][1]*spin_m[1]; 
	printf( "| sigmaDOTm> %e %e \n", v1[0].real, v1[1].real); 


	p = probability(v1, v1); 	printf( "|<v1|v1>|^2  %e \n", p); 
	p = probability(v2, v2); 	printf( "|<v2|v2>|^2  %e \n", p); 
	p = probability(v1, v2); 	printf( "|<v1|v2>|^2  %e \n", p); 

	p = m_x*n_x + m_y*n_y + m_z*n_z; 	printf( "%e n*m \n",  p/2+0.5); 
	p = (1.0+cos(alpha)) /2.0; 	printf( "%e  cos \n", p); 
} 




int main()//gravity curves
{
	init_system();
	
	int h2=1;
	double z=0.5,dy=0,dx=1e-2,dy2=0,dy3=0,dy4=0,dx2=0;
	for(int x2=0;x2<1000;x2+=h2)
	{
		double x=1+x2*dx;
		
		//F=Dx Hooks law
		double du= 1e1/(3.0*x*x);
		dy+=sqrt(du*du-dx*dx);
		int y=dy;	if(y<0) y=0;	if(y>400) y=400;
		pixel(x2,500-y,0xff0000);

		du= 1e1/(3.0*x*x*x);
		dy3+=sqrt(du*du-dx*dx);
		y=dy3;	if(y<0) y=0;	if(y>400) y=400;
		pixel(x2,500-y,0x0000ff);


		double rs=1e0;
		dy4=2*sqrt(rs*((x)-rs))*5e1;//flamm
		y=dy4;	if(y<0) y=0;	if(y>400) y=400;
		pixel(x2,500-y,0x00ffff);

		rs=1e-2;
		x=rs*2.00001+x2*1e-5;
		dx2+=sqrt(1.0-rs/(x-rs));
		dy2+=(1.0/sqrt(1.0-rs/(x-rs)));
		y=dy2*1e-1;	if(y<0) y=0;	if(y>400) y=400;
		pixel(dx2,500-y,0x00ff00);
//printf("%e \n",1.0/sqrt(1.0-rs/(x)));		
//printf("%e %e\n",dx2,x);		
//printf("%e \n",dy2);		
	}

	XFlush(dpy);
	getchar();

	return 0;
}


#endif


