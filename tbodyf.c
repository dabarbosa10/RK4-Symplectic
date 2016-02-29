#include <stdio.h>
#include <math.h>

void leapfrog_step(double delta_t, double t, double *q1, double *p1,double *q2, double*p2);
void leapfrog_step1(double delta_t, double t, double *q1, double *p1, double *q2, double *p2);
void RK4_step(double delta_t, double t, double *q1, double *p1, double *q2, double *p2);
double deriv_q1(double t, double q1, double p1, double q2, double p2);
double deriv_p1(double t, double q1, double p1, double q2, double p2);
double deriv_q2(double t, double q1, double p1, double q2, double p2);
double deriv_p2(double t, double q1, double p1, double q2, double p2);

int main(){

   double q1_RK4, p1_RK4, q2_RK4, p2_RK4;
  double  q1_LF;
  double  p1_LF;
  double  q2_LF;
  double  p2_LF;

  double t;

  double T=2000.0;
  double delta_t=0.006;
  int n_step = (int)(T/delta_t);
  int i,j;
  
  for(j=0;j<=80;j++){
  t=0.0;
  q1_RK4=0.3555;
  p1_RK4=0.0;
  q2_RK4=0.0005+0.045*j;
  p2_RK4=0.0;
  q1_LF=0.3555;
  p1_LF=0.0;
  q2_LF=0.0005+0.045*j;
  p2_LF=0.0;
  for(i=0;i<n_step;i++){    
     if( p1_LF>=-0.0009 && p1_LF<=0.0009 && p1_RK4>=-0.0009 && p1_RK4<=0.0009 ){
    printf("  %f \t %f \t  %f \t %f \n", q2_LF, p2_LF,q2_RK4, p2_RK4);}
    RK4_step(delta_t, t, &q1_RK4, &p1_RK4, &q2_RK4, &p2_RK4);
    leapfrog_step(delta_t, t, &q1_LF, &p1_LF, &q2_LF, &p2_LF);
    leapfrog_step1(delta_t, t, &q1_LF, &p1_LF, &q2_LF, &p2_LF);
    t += delta_t;
  }
 
}

}

double deriv_q1(double t, double q1, double p1,double q2, double p2 ){
  return p1;
}

double deriv_p1(double t, double q1, double p1, double q2, double p2){
double r1=(4*q1*q1)+1;
double pot=3/2;
return (  -(2*q1)/ pow(r1, 1.5)  );

}

double deriv_q2(double t, double q1, double p1, double q2, double p2){
return (p2); 
}

double deriv_p2(double t, double q1, double p1, double q2, double p2){
double r2=q1-q2;
double r3=q1+q2;
double pot1=1.5;
double fu1= r2*r2+0.25;
double fu2= r3*r3+0.25;
return ( (r2/pow(fu1,1.5)) - (r3/pow(fu2,1.5))                  );
}

void leapfrog_step(double delta_t, double t, double *q1, double *p1, double *q2, double *p2){
   double q1_in;
  double p1_in;
  double q2_in;
  double p2_in;
  q1_in = *q1;
  p1_in = *p1;
  q2_in = *q2;
  p2_in = *p2;

  /*kick*/
  p1_in += 0.5 * deriv_p1(t,q1_in,p1_in, q2_in, p2_in) * delta_t;
  /*drift*/
  q1_in += 1.0 * p1_in * delta_t;
  /*kick*/
  p1_in += 0.5 * deriv_p1(t, q1_in, p1_in, q2_in, p2_in) * delta_t;


  *q1 = q1_in;
  *p1= p1_in;
}

void leapfrog_step1(double delta_t, double t, double *q1, double *p1, double *q2, double *p2){
  double q1_in;
  double p1_in;
  double q2_in;
  double p2_in;
  q1_in = *q1;
  p1_in = *p1;
  q2_in = *q2;
  p2_in = *p2;

  /*kick*/
  p2_in += 0.5 * deriv_p2(t,q1_in,p1_in, q2_in, p2_in) * delta_t;
  /*drift*/
  q2_in += 1.0 * p2_in * delta_t;
  /*kick*/
  p2_in += 0.5 * deriv_p2(t, q1_in, p1_in,q2_in,p2_in) * delta_t;


  *q2 = q2_in;
  *p2= p2_in;
}

void RK4_step(double delta_t, double t, double *q1, double *p1, double *q2, double *p2){
  double k1, k2, k3, k4;
  double l1, l2, l3, l4;
  double m1, m2, m3, m4;
  double n1, n2, n3, n4;
  double q1_in;
  double p1_in;
  double q2_in;
  double p2_in;
  q1_in = *q1;
  p1_in = *p1;
  q2_in = *q2;
  p2_in = *p2;



  k1 = deriv_q1(t,q1_in,p1_in,q2_in,p2_in);  
  l1 = deriv_q2(t,q1_in,p1_in,q2_in,p2_in);  
  m1 = deriv_p1(t,q1_in,p1_in,q2_in,p2_in);  
  n1 = deriv_p2(t,q1_in,p1_in,q2_in,p2_in);  

  k2 = deriv_q1(t + delta_t*0.5, q1_in + k1*delta_t*0.5, p1_in + m1*delta_t*0.5, q2_in + l1*delta_t*0.5, p2_in +n1*delta_t*0.5);
  l2 = deriv_q2(t + delta_t*0.5, q1_in + k1*delta_t*0.5, p1_in + m1*delta_t*0.5, q2_in + l1*delta_t*0.5, p2_in +n1*delta_t*0.5);
  m2 = deriv_p1(t + delta_t*0.5, q1_in + k1*delta_t*0.5, p1_in + m1*delta_t*0.5, q2_in + l1*delta_t*0.5, p2_in +n1*delta_t*0.5);
  n2 = deriv_p2(t + delta_t*0.5, q1_in + k1*delta_t*0.5, p1_in + m1*delta_t*0.5, q2_in + l1*delta_t*0.5, p2_in +n1*delta_t*0.5);

  k3 = deriv_q1(t + delta_t*0.5, q1_in + k2*delta_t*0.5, p1_in + m2*delta_t*0.5, q2_in + l2*delta_t*0.5, p2_in + n2*delta_t*0.5);
  l3 = deriv_q2(t + delta_t*0.5, q1_in + k2*delta_t*0.5, p1_in + m2*delta_t*0.5, q2_in + l2*delta_t*0.5, p2_in + n2*delta_t*0.5);
  m3 = deriv_p1(t + delta_t*0.5, q1_in + k2*delta_t*0.5, p1_in + m2*delta_t*0.5, q2_in + l2*delta_t*0.5, p2_in + n2*delta_t*0.5);
  n3 = deriv_p2(t + delta_t*0.5, q1_in + k2*delta_t*0.5, p1_in + m2*delta_t*0.5, q2_in + l2*delta_t*0.5, p2_in + n2*delta_t*0.5);

  k4 = deriv_q1(t + delta_t, q1_in + k3*delta_t, p1_in + m3*delta_t, q2_in + l3*delta_t, p2_in + n3*delta_t);
  l4 = deriv_q2(t + delta_t, q1_in + k3*delta_t, p1_in + m3*delta_t, q2_in + l3*delta_t, p2_in + n3*delta_t);
  m4 = deriv_p1(t + delta_t, q1_in + k3*delta_t, p1_in + m3*delta_t, q2_in + l3*delta_t, p2_in + n3*delta_t);
  n4 = deriv_p2(t + delta_t, q1_in + k3*delta_t, p1_in + m3*delta_t, q2_in + l3*delta_t, p2_in + n3*delta_t);

  q1_in += (k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0)*delta_t;
  q2_in += (l1/6.0 + l2/3.0 + l3/3.0 + l4/6.0)*delta_t;
  p1_in += (m1/6.0 + m2/3.0 + m3/3.0 + m4/6.0)*delta_t;
  p2_in += (n1/6.0 + n2/3.0 + n3/3.0 + n4/6.0)*delta_t;

  *q1 = q1_in;
  *q2 = q2_in;
  *p1 = p1_in;
  *p2 = p2_in;
}

