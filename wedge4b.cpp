
//g++ x.cpp  -O3 


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


float delta_e=0.1;

//#define float1 float
#define float1 double

#define iradian (M_PI/180.0)
#define frnd(n) (frand()*(float)n)





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

vec3 operator - () {x=-x;y=-y;z=-z;return *this;};
};

float1 dot(vec3 v1,vec3 v2) { return (v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);}
float1 length(vec3 v1) 
{
	float1 t1=dot(v1,v1);
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


#if 0
mat*mat
mat3[a][b]=mat1[a][i]*mat2[i][b]=
mat3[a][b]=mat1[a][1]*mat2[1][b] + mat1[a][2]*mat2[2][b] + mat1[a][3]*mat2[3][b] + mat1[a][4]*mat2[4][b] 

	mo._11=ma_11*mb_11 + ma_12*mb_21 + ma_13*mb_31 + ma_14*mb_41 ;
	mo._12=ma_11*mb_12 + ma_12*mb_22 + ma_13*mb_32 + ma_14*mb_42 ;
	mo._13=ma_11*mb_13 + ma_12*mb_23 + ma_13*mb_33 + ma_14*mb_43 ;
	mo._14=ma_11*mb_14 + ma_12*mb_24 + ma_13*mb_34 + ma_14*mb_44 ;

	mo._21=ma_21*mb_11 + ma_22*mb_21 + ma_23*mb_31 + ma_24*mb_41 ;
	mo._22=ma_21*mb_12 + ma_22*mb_22 + ma_23*mb_32 + ma_24*mb_42 ;
	mo._23=ma_21*mb_13 + ma_22*mb_23 + ma_23*mb_33 + ma_24*mb_43 ;
	mo._24=ma_21*mb_14 + ma_22*mb_24 + ma_23*mb_34 + ma_24*mb_44 ;

	mo._31=ma_31*mb_11 + ma_32*mb_21 + ma_33*mb_31 + ma_34*mb_41 ;
	mo._32=ma_31*mb_12 + ma_32*mb_22 + ma_33*mb_32 + ma_34*mb_42 ;
	mo._33=ma_31*mb_13 + ma_32*mb_23 + ma_33*mb_33 + ma_34*mb_43 ;
	mo._34=ma_31*mb_14 + ma_32*mb_24 + ma_33*mb_34 + ma_34*mb_44 ;

	mo._41=ma_41*mb_11 + ma_42*mb_21 + ma_43*mb_31 + ma_44*mb_41 ;
	mo._42=ma_41*mb_12 + ma_42*mb_22 + ma_43*mb_32 + ma_44*mb_42 ;
	mo._43=ma_41*mb_13 + ma_42*mb_23 + ma_43*mb_33 + ma_44*mb_43 ;
	mo._44=ma_41*mb_14 + ma_42*mb_24 + ma_43*mb_34 + ma_44*mb_44 ;
#endif


//##############################################################

struct tensor_t
{
	float1 xx;
	float1 xy;
	float1 xz;
	float1 yx;
	float1 yy;
	float1 yz;
	float1 zx;
	float1 zy;
	float1 zz;
	
	tensor_t()
	{
		xx=0.0;
		xy=0.0;
		xz=0.0;
		yx=0.0;
		yy=0.0;
		yz=0.0;
		zx=0.0;
		zy=0.0;
		zz=0.0;
	}
	
	void scale(float1 scl)
	{
		xx*=scl;
		xy*=scl;
		xz*=scl;
		yx*=scl;
		yy*=scl;
		yz*=scl;
		zx*=scl;
		zy*=scl;
		zz*=scl;
	}
	tensor_t operator +=(tensor_t v)
	{
		xx+=v.xx;
		xy+=v.xy;
		xz+=v.xz;
		yx+=v.yx;
		yy+=v.yy;
		yz+=v.yz;
		zx+=v.zx;
		zy+=v.zy;
		zz+=v.zz;

		return *this;
	}
	tensor_t operator -=(tensor_t v)
	{
		xx-=v.xx;
		xy-=v.xy;
		xz-=v.xz;
		yx-=v.yx;
		yy-=v.yy;
		yz-=v.yz;
		zx-=v.zx;
		zy-=v.zy;
		zz-=v.zz;
		
		return *this;
	}
//https://en.wikipedia.org/wiki/Rule_of_Sarrus
//xx xy xz    a,b,c  aef  -bdf   cde
//yx yy yz    d,e,f   hi    gi    gh
//zx zy zz    g,h,i
	float determinant()
	{
/*		printf("dt %e %e %e %e %e %e \n", // OKAY
			 	xx*yy*zz ,-xx*yz*zy ,
				xy*yz*zx ,-xy*yx*zz , 
				xz*yx*zy ,-xz*yy*zx);*/
				
		return 	xx*yy*zz -xx*yz*zy +    // (aei - afh)
				xy*yz*zx -xy*yx*zz +    //-(bdi - bfg) = (bfg - bdi)
				xz*yx*zy -xz*yy*zx;     // (cdh - ceg)
	}
};


struct bivector
{
    tensor_t tensor;

    bivector()
    {
    }
    
    //https://en.wikipedia.org/wiki/Exterior_algebra
    //The exterior product is sometimes called the outer product, although this can also refer to the tensor product of vectors.
    static bivector geometric_product(vec3 a, vec3 b)   
    {
	    bivector e;
	    
	    e.tensor=outer_product(a,b);
	    return e;
    }
    
//   https://en.wikipedia.org/wiki/Outer_product
	static tensor_t outer_product(vec3 a, vec3 b)
	{
		tensor_t e;
		e.xx=a.x*b.x;
		e.xy=a.x*b.y;
		e.xz=a.x*b.z;

		e.yx=a.y*b.x;
		e.yy=a.y*b.y;
		e.yz=a.y*b.z;// u2v3

		e.zx=a.z*b.x;
		e.zy=a.z*b.y;// u3v2
		e.zz=a.z*b.z;
// u V v = (u2v3 - u3v2)(e2 V e3) + (u3v1 - u1v3)(e3 V e1) + (u1v2 - u2v1)(e1 V e2) 
// x=tensor.yz-tensor.zy  OK
//tehat lehet wedge egy vektor(bivector bazissal)  ha a determinansokbol epitjuk fel
		return e;
	}

//Hodge star operators
//xx xy xz  1  0  0 
//yx yy yz  0  1  0 
//zx zy zz  0  0  1 
//https://en.wikipedia.org/wiki/Exterior_algebra#Inner_product
    operator float1()  { return tensor.xx + tensor.yy + tensor.zz; } 
    
//xx xy xz  0  1 -1 
//yx yy yz -1  0  1 
//zx zy zz  1 -1  0 
//https://en.wikipedia.org/wiki/Exterior_algebra#Cross_and_triple_products    
    operator vec3()   { return vec3(tensor.yz-tensor.zy, tensor.zx-tensor.xz, tensor.xy-tensor.yx); } 
/*    operator vec3()   { return 	vec3(
    					tensor.xx + tensor.yz - tensor.zy, // a cross kijon a (ab-ba)bol, ha ferden van a dot() REF5
    				   -tensor.xz + tensor.yy + tensor.zx, //NA nem jo igy, mert a dot belelog!!
    					tensor.xy - tensor.yx + tensor.zz); } */


	float1 determinant()
	{
		return tensor.determinant();
	}

    bivector operator *( float1 scale_value)
    {
	    bivector e=*this;
        e.tensor.scale(scale_value);

	    return e;
	}
    bivector operator +( bivector parm)
    {
	    bivector e=*this;
        e.tensor += parm.tensor;

        return e;
    }
    bivector operator -(bivector parm)
    {
	    bivector e=*this;
        e.tensor -= parm.tensor;

        return e;
    }
    void dump(const char *ss)
    {
    	printf("%s \n",ss);
	    printf("%e %e %e \n",tensor.xx,tensor.xy,tensor.xz);
	    printf("%e %e %e \n",tensor.yx,tensor.yy,tensor.yz);
	    printf("%e %e %e \n\n",tensor.zx,tensor.zy,tensor.zz);

/*		float1 vx=tensor.xx + tensor.xy + tensor.xz;  // NA        a cross kijon a (ab-ba)bol, ha ferden van a dot() REF5
		float1 vy=tensor.yx + tensor.yy + tensor.yz;
		float1 vz=tensor.zx + tensor.zy + tensor.zz;
	    printf("%e, %e ,%e ,\n",vx,vy,vz);

		vx=tensor.xx + tensor.yx + tensor.zx;
		vy=tensor.xy + tensor.yy + tensor.zy;
		vz=tensor.xz + tensor.yz + tensor.zz;
	    printf("%e, %e ,%e ,\n\n",vx,vy,vz);*/
	}
};


int main()
{
#if 0
	for(int a=0;a<4;a++)// Fab
	{
		for(int b=0;b<4;b++) 
		{
			printf("F[%d][%d]=dA%d.%d - dA%d.%d \n",a,b,a,b,b,a);
		}
		printf("\n");
	}	
#endif

	
#if 1
    vec3 a = vec3(5,23,2);
    vec3 b = vec3(12,-5,-15);

    bivector ab=bivector::geometric_product(a, b);
    bivector ba=bivector::geometric_product(b, a);//  transpose matrix

    bivector abPLUSba = (ab + ba)*0.5;
    bivector abMINUSba = (ab - ba)*0.5; // wedge !!!
    bivector ab2 = abPLUSba + abMINUSba;

ab.dump("ab");
ab2.dump("ba2 (ab+ba + ab-ba)");
ba.dump("ba");                   //  transpose matrix
abPLUSba.dump("ab+ba");
abMINUSba.dump("ab-ba");

    float1 hodge1=(float1)abPLUSba;
    vec3 hodge1v=(vec3)abPLUSba;
	printf(" abPLUSba scalar %e vector %e %e %e \n",hodge1,hodge1v.x,hodge1v.y,hodge1v.z);
    float1 hodge2=(float1)abMINUSba;
    vec3 hodge2v=(vec3)abMINUSba;
	printf(" abMINUSba scalar %e vector %e %e %e \n",hodge2,hodge2v.x,hodge2v.y,hodge2v.z);
    
    
	printf("The interior product \n");
	float1 dot3=(dot(a+b,a+b)-dot(a,a)- dot(b,b))*0.5;
	printf("  %e == %e \n",hodge1,dot3);
	printf("  %e == %e == %e\n",hodge1,(float1)ab,(float1)ba);

	printf("The exterior product \n");
	printf("  %e %e %e != \n  %e %e %e \n",
		((vec3)ab).x,((vec3)ab).y,((vec3)ab).z,
		((vec3)ba).x,((vec3)ba).y,((vec3)ba).z);

	vec3 nor=cross(a,b);
	printf("The cross product \n");
	printf("  %e %e %e \n",		nor.x,nor.y,nor.z);

//ab.tensor.zx+=15;
	printf("det %e \n",ab.determinant());
	printf("    %e \n",ba.determinant());
	printf("    %e \n",abPLUSba.determinant());
	printf("    %e \n",abPLUSba.determinant());


    float1 check_1 = (float1)ab - (float1)abPLUSba;         // 0
    vec3 check_2 = (vec3)ab - (vec3)abMINUSba;     // (0,0,0)

	printf("%e \n",check_1);
	printf("%e %e %e \n",check_2.x,check_2.y,check_2.z);
#endif	

return 0;
}

/*
https://en.wikipedia.org/wiki/Bivector#The_interior_product
https://en.wikipedia.org/wiki/Geometric_algebra
https://en.wikipedia.org/wiki/Exterior_algebra



ab=(ab+ba)/2 + (ab-ba)/2
ab=dot(a,b) + (a cross b)
ab
(a1)(b1+b2+b3)
(a2)
(a3)
a1b1 , a1b2 , a1b3
a2b1 , a2b2 , a2b3
a3b1 , a3b2 , a3b3

(ab+ba)/2
a1b1 + b1a1 , a1b2 + b1a2 , a1b3 + b1a3  /2
a2b1 + b2a1 , a2b2 + b2a2 , a2b3 + b2a3  /2
a3b1 + b3a1 , a3b2 + b3a2 , a3b3 + b3a3  /2

(ab-ba)/2 
a1b1 - b1a1 , a1b2 - b1a2 , a1b3 - b1a3  /2
a2b1 - b2a1 , a2b2 - b2a2 , a2b3 - b2a3  /2
a3b1 - b3a1 , a3b2 - b3a2 , a3b3 - b3a3  /2

(ab-ba)/2 
          0 , a1b2 - b1a2 , a1b3 - b1a3  /2
a2b1 - b2a1 ,           0 , a2b3 - b2a3  /2
a3b1 - b3a1 , a3b2 - b3a2 ,           0  /2


dot product    
xx xy xz  1  0  0 
yx yy yz  0  1  0 
zx zy zz  0  0  1 
    
cross product    
xx xy xz  0  z -y 
yx yy yz -z  0  x 
zx zy zz  y -x  0 



(ab+ba)/2
dot
(2*a1b1 + 2*a2b2 + 2*a3b3 )/2 = dot(a,b)

cross
x = ((a2b3 + b2a3) - (a3b2 + b3a2)) /2 =0
y = ((a3b1 + b3a1) - (a1b3 + b1a3)) /2 =0
z = ((a1b2 + b1a2) - (a2b1 + b2a1)) /2 =0

(ab-ba)/2
dot
(a1b1 - b1a1 + a2b2 - b2a2 + a3b3 - b3a3 )/2 = 0

cross
x = ((a2b3 - b2a3) - (a3b2 - b3a2)) /2
y = ((a3b1 - b3a1) - (a1b3 - b1a3)) /2
z = ((a1b2 - b1a2) - (a2b1 - b2a1)) /2

x = (a2b3 - b2a3 - a3b2 + b3a2) /2
y = (a3b1 - b3a1 - a1b3 + b1a3) /2
z = (a1b2 - b1a2 - a2b1 + b2a1) /2

x = (a2b3 - b2a3 - a3b2 + b3a2) /2
y = (a3b1 - b3a1 - a1b3 + b1a3) /2
z = (a1b2 - b1a2 - a2b1 + b2a1) /2

x = (2*a2b3 - 2*a3b2 ) /2
y = (2*a3b1 - 2*a1b3 ) /2
z = (2*a1b2 - 2*a2b1 ) /2

x = a2b3 - a3b2
y = a3b1 - a1b3
z = a1b2 - a2b1


hogyan lehet termeszetessebbe tenni ezt a keveredest, 
hogyan kapjuk meg a keveredest automatikusan?
ab=(ab+ba)/2 + (ab-ba)/2
ab=dot(a,b) + (a wedge b)

dual base ssolution, grassmann algebra
e1 {1,0,0}
e2 {0,1,0}
e3 {0,0,1}

a = a1e1 + a2e2 + a3e3
b = b1e1 + b2e2 + b3e3

ab = (a1e1 + a2e2 + a3e3)*(b1e1 + b2e2 + b3e3)
ab = 
a1e1*b1e1 + a2e2*b1e1 + a3e3*b1e1 +
a1e1*b2e2 + a2e2*b2e2 + a3e3*b2e2 +
a1e1*b3e3 + a2e2*b3e3 + a3e3*b3e3 

ab = 
a1b1*e1e1 + a2b1*e2e1 + a3b1*e3e1 +    e1e1=1
a1b2*e1e2 + a2b2*e2e2 + a3b2*e3e2 +
a1b3*e1e3 + a2b3*e2e3 + a3b3*e3e3 

ab = 
a1b1      + a2b1*e2e1 + a3b1*e3e1 +
a1b2*e1e2 + a2b2      + a3b2*e3e2 +
a1b3*e1e3 + a2b3*e2e3 + a3b3     

ab = 
a1b1 + a2b2 + a3b3   +
a2b3*e2e3 + a3b2*e3e2 +
a3b1*e3e1 + a1b3*e1e3 +
a1b2*e1e2 + a2b1*e2e1 +

a sorrend nem szamit, mivel az e2e3 mondja meg, hogy ez pl X
mivel e1e2 = -e2e1 a cross szerint

ab = 
a1b1 + a2b2 + a3b3   +
(a2b3 - a3b2)*e2e3 +
(a3b1 - a1b3)*e1e3 +
(a1b2 - a2b1)*e1e2 +

ab = dot(a,b) + (a wedge b)
                (2blade)



hodge *
c = -Ai
i=e1e2e3

Ai=
(a2b3 - a3b2)*e2e3 *e1e2e3+       e1^2=e1e1=e2e2=e3e3=1   
(a3b1 - a1b3)*e1e3 *e1e2e3+
(a1b2 - a2b1)*e1e2 *e1e2e3+

(a2b3 - a3b2)*e1(e2e3)^2+      (e1e2)^2=-1    mivel a sorrend valtozik,(cross prduct)
(a3b1 - a1b3)*e2(e1e3)^2+
(a1b2 - a2b1)*e3(e1e2)^2+

-(a2b3 - a3b2)*e1+ 
-(a3b1 - a1b3)*e2+
-(a1b2 - a2b1)*e3+
x=-(a2b3 - a3b2)
y=-(a3b1 - a1b3)
z=-(a1b2 - a2b1)
=c = -(*A)
 
--------------------------------------------

references
cross(a,b)
a2b3 - a3b2  
a3b1 - a1b3
a1b2 - a2b1


*/




