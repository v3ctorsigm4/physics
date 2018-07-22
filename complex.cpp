

#if 1

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
void add_polarizer3(vec2 *v,vec2 *residuo,float phase)
{
	vec2 pol;

	add_amp(&pol,phase,1.0);
	float amp=dot(*v,pol);      //v'=pol><pol|v>
    v->x = pol.x*amp;    
    v->y = pol.y*amp;
    residuo->x=pol.x - v->x;
    residuo->y=pol.y - v->y;
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
void add_polarizer3(complex *v,complex *residuo,float phase)
{
	complex pol,pol2;
	add_amp(&pol,phase,1.0);
	pol2=pol;
	pol2.conjugate();
	complex amp=(*v)*pol2; amp.img=0;      //v'=pol><pol|v>
    *v = pol*amp;    
    residuo->x=pol.x - v->x;
    residuo->y=pol.y - v->y;
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

;
