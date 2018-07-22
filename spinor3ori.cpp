
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




