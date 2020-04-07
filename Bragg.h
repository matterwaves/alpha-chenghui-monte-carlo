#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_spline.h"


int m0, HSize, type, n, picture;
double Omega1, Omega2, Delta, tc, kc, phase, omega_m, phic, delta, rampon, rampoff, ti, pulsewidth;
gsl_interp_accel *acc, *acc1;
gsl_spline *spline, *spline1;

//Complex number struct definition

struct complex
{
  double real,imag;
};

struct complex I =
{
  .real = 0,
  .imag = 1
};

struct complex complex_add(struct complex a, struct complex b)
{
  struct complex ans;
  ans.real = a.real+b.real;
  ans.imag = a.imag+b.imag;
  return ans;
}

struct complex complex_mul_real(struct complex a, double c)
{
  struct complex ans;
  ans.real = c*a.real;
  ans.imag = c*a.imag;
  return ans;
}

struct complex complex_mul(struct complex a, struct complex b)
{
  struct complex ans;
  ans.real = a.real*b.real-a.imag*b.imag;
  ans.imag = a.real*b.imag+a.imag*b.real;
  return ans;
}

struct complex complex_exp(double x)
{
  struct complex ans;
  ans.real = cos(x);
  ans.imag = sin(x);
  return ans;
}

//Complex array operation
void complex_vector_add(struct complex *a, const struct complex *b, int size)
{
  int i;
  for (i = 0; i < size; i++)
    a[i] = complex_add(a[i],b[i]);
}

void complex_vector_mul_real(struct complex *a, double c, int size)
{
  int i;
  for (i = 0; i < size; i++)
    a[i] = complex_mul_real(a[i],c);
}

void complex_vector_mul(struct complex *a, struct complex b, int size)
{
  int i;
  for (i = 0; i < size; i++)
    a[i] = complex_mul(a[i],b);
}

void complex_vector_copy(struct complex *copy, const struct complex *a, int size)
{
  int i;
  for (i = 0; i < size; i++)
    copy[i] = a[i];
}

//Science section I (Bragg diffraction)

double pulse(double t)
{
  if (type > 1)
  {
      if (t < ti+rampon)
          return (t-ti)/rampon;
      if (t > ti+pulsewidth-rampoff)
          return -(t-ti-pulsewidth)/rampoff;
      return 1;
  }
  return gsl_spline_eval(spline, t, acc);
}

double ramprate(double t)
{
  if (t < ti+rampon)
    return 0;
  if (t > ti+pulsewidth-rampoff)
    return 0;
  return gsl_spline_eval(spline1, t, acc1);
}

void Hamiltonian(struct complex * dpsi, const struct complex * psi, double t)
{
  int i;
  double alpha1, alpha2, beta;

  m0 = -n/2-n-6;
  if (picture == 0) //Schrodinger picture
  {
    if (type == 0 || type == 2) //single frequency
    {
      //alpha1 = (Omega1*Omega1+Omega2*Omega2)/4.0/Delta;
      alpha1 = 0;
      beta = Omega1*Omega2/4.0/Delta;
      for (i = 0; i < HSize; i++)
      {
	dpsi[i] = complex_mul_real(psi[i],(kc+(m0+i)*2)*(kc+(m0+i)*2)+alpha1*pulse(t));
	if (i > 0)
	  dpsi[i] = complex_add(dpsi[i],complex_mul(complex_mul_real(complex_exp(phic-delta*t),beta*pulse(t)),psi[i-1]));
	if (i < HSize-1)
	  dpsi[i] = complex_add(dpsi[i],complex_mul(complex_mul_real(complex_exp(-phic+delta*t),beta*pulse(t)),psi[i+1]));
      }
    }
    else //multi frequencies
    {
      //alpha1 = (2*Omega1*Omega1+Omega2*Omega2)/4.0/Delta;
      //alpha2 = Omega1*Omega1/2.0/Delta;
      alpha1 = 0;
      alpha2 = 0;
      beta = Omega1*Omega2/4.0/Delta;
      for (i = 0; i < HSize; i++)
      {
	dpsi[i] = complex_mul_real(psi[i],(kc+(m0+i)*2)*(kc+(m0+i)*2)+alpha1*pulse(t)+alpha2*cos(2*omega_m*t)*pulse(t));
	if (i > 0)
	  dpsi[i] = complex_add(dpsi[i],complex_mul(complex_mul_real(complex_exp(phic-delta*t),2*beta*pulse(t)*cos(omega_m*t)),psi[i-1]));
	if (i < HSize-1)
	  dpsi[i] = complex_add(dpsi[i],complex_mul(complex_mul_real(complex_exp(-phic+delta*t),2*beta*pulse(t)*cos(omega_m*t)),psi[i+1]));
      }
    }
  }
  else //Interaction picture
  {
    if (type == 0 || type == 2) //single frequency
    {
      //alpha1 = (Omega1*Omega1+Omega2*Omega2)/4.0/Delta;
      alpha1 = 0;
      beta = Omega1*Omega2/4.0/Delta;
      for (i = 0; i < HSize; i++)
      {
        dpsi[i].real = 0;
        dpsi[i].imag = 0;
        if (i > 0)
          dpsi[i] = complex_add(dpsi[i],complex_mul(complex_mul_real(complex_exp(phic-delta*t+4*(kc+2*(m0+i)-1)*t),beta*pulse(t)),psi[i-1]));
        if (i < HSize-1)
          dpsi[i] = complex_add(dpsi[i],complex_mul(complex_mul_real(complex_exp(-phic+delta*t-4*(kc+2*(m0+i)+1)*t),beta*pulse(t)),psi[i+1]));
      }
    }
    else //multi frequencies
    {
      //alpha1 = (2*Omega1*Omega1+Omega2*Omega2)/4.0/Delta;
      //alpha2 = Omega1*Omega1/2.0/Delta;
      alpha1 = 0;
      alpha2 = 0;
      beta = Omega1*Omega2/4.0/Delta;
      for (i = 0; i < HSize; i++)
      {
        dpsi[i].real = 0;
        dpsi[i].imag = 0;
        if (i > 0)
          dpsi[i] = complex_add(dpsi[i],complex_mul(complex_mul_real(complex_exp(phic-delta*t+4*(kc+2*(m0+i)-1)*t),2*beta*pulse(t)*cos(omega_m*t)),psi[i-1]));
        if (i < HSize-1)
          dpsi[i] = complex_add(dpsi[i],complex_mul(complex_mul_real(complex_exp(-phic+delta*t-4*(kc+2*(m0+i)+1)*t),2*beta*pulse(t)*cos(omega_m*t)),psi[i+1]));
      }
    }
  }
  complex_vector_mul(dpsi, complex_mul_real(I,-1), HSize);
}

//Scientific computing section (Runge-Kutta)
void Runge_Kutta(struct complex * psi_out, double t0, double * aoi, int aoi_size, double stepsize)
{
  double t, amp, dt;
  struct complex *k1, *k2, *k3, *k4;
  struct complex *tmp, *psi;
  int i, j;

  t = t0;
  i = 1;
  psi = (struct complex *)malloc(sizeof(struct complex)*HSize);
  k1 = (struct complex *)malloc(sizeof(struct complex)*HSize);
  k2 = (struct complex *)malloc(sizeof(struct complex)*HSize);
  k3 = (struct complex *)malloc(sizeof(struct complex)*HSize);
  k4 = (struct complex *)malloc(sizeof(struct complex)*HSize);
  tmp = (struct complex *)malloc(sizeof(struct complex)*HSize);
  for (j = 0; j < HSize; j++)
    psi[j] = psi_out[j];

  while (i < aoi_size)
  {
    if (t+stepsize > aoi[i]+t0)
      dt = aoi[i]+t0-t;
    else
      dt = stepsize;
    /*
    if (i == 5000)
    {
      cout.precision(12);
      for (j = 0; j < HSize; j++)
      {
	cout << norm(H(j,0));
	for (k = 1; k < HSize; k++) cout << '\t' << norm(H(j,k));
	cout << endl;
      }
      cout << endl;
      for (j = 0; j < HSize; j++)
      {
	cout << arg(H(j,0));
	for (k = 1; k < HSize; k++) cout << '\t' << arg(H(j,k));
	cout << endl;
      }
      cout << endl;
    }
    */
    //k1 = -iH(t, psi)
    Hamiltonian(k1, psi, t);
    //printf("k1=%f\n",k1[9].imag);

    //k2 = -iH(t+dt/2, psi+k1*dt/2);
    complex_vector_copy(tmp, k1, HSize);
    complex_vector_mul_real(tmp, dt/2, HSize);
    complex_vector_add(tmp, psi, HSize);
    Hamiltonian(k2, tmp, t+dt/2);
    //printf("k2=%f\n",k2[9].imag);

    //k3 = -iH(t+dt/2, psi+k2*dt/2);
    complex_vector_copy(tmp, k2, HSize);
    complex_vector_mul_real(tmp, dt/2, HSize);
    complex_vector_add(tmp, psi, HSize);
    Hamiltonian(k3, tmp, t+dt/2);
    //printf("k3=%f\n",k3[9].imag);

    //k4 = -iH(t+dt, psi+k3*dt);
    complex_vector_copy(tmp, k3, HSize);
    complex_vector_mul_real(tmp, dt, HSize);
    complex_vector_add(tmp, psi, HSize);
    Hamiltonian(k4, tmp, t+dt);
    //printf("k4=%f\n",k4[9].imag);

    //psi(n+1) = psi(n)+dt/6*(k1+2*k2+2*k3+k4);
    complex_vector_mul_real(k2, 2, HSize);
    complex_vector_mul_real(k3, 2, HSize);
    complex_vector_copy(tmp, k4, HSize);
    complex_vector_add(tmp, k3, HSize);
    complex_vector_add(tmp, k2, HSize);
    complex_vector_add(tmp, k1, HSize);
    complex_vector_mul_real(tmp, dt/6, HSize);
    complex_vector_add(psi, tmp, HSize);

    if (type > 1)
    {
        kc += ramprate(t)/4*dt;
    }
    //printf("t=%f, delta=%f, kc=%f\n", t, delta, kc);
    t = t+dt;
    /*
    amp = 0;

    for (j = 0; j < HSize; j++)
    {
      amp += psi[j].real*psi[j].real+psi[j].imag*psi[j].imag;
      printf("%f, %d, %f, %f\n",t,j,psi[j].real,psi[j].imag);
    }
    printf("%f, %f\n", t, amp);

    printf("%f, %f\n", t, amp);
    if (i%100 == 0)
      printf("%f, %f+%fi\n", t, psi[10].real, psi[10].imag);
    */
    if (dt < stepsize)
    {
      for (j = 0; j < HSize; j++)
	psi_out[i*HSize+j] = psi[j];
      i = i+1;
    }
  }
  free(k1);
  free(k2);
  free(k3);
  free(k4);
  free(tmp);
}

//Main Functions
void Bragg(double * psi_real, double * psi_imag, double init_v, int picture0, int type0, int n0, double Omega10, double Omega20, double Delta0, double delta0, double pulsewidth0, double rampon0, double rampoff0, double omega_m0, double t0, double phic0, double stepsize, double * pulseshape, double * rampshape, int pulseshape_size, int rampshape_size, double * aoi, int aoi_size)
{
  struct complex * psi;
  int i, index;
  double * tlist;
  //Get paramters
  ti = t0;
  pulsewidth = pulsewidth0;
  Omega1 = Omega10;
  Omega2 = Omega20;
  Delta = Delta0;
  delta = delta0;
  omega_m = omega_m0;
  phic = phic0;
  picture = picture0;
  type = type0;
  rampon = rampon0;
  rampoff = rampoff0;
  n = n0;
  HSize = n*3+11;
  m0 = -n/2-n-6;
  kc = init_v-floor(init_v/2)*2;
  index = floor(init_v/2)-m0;

  //pulsewidth0 = 1.2;
  //Pulse Shape Interpolation
  tlist = (double *)malloc(sizeof(double)*pulseshape_size);
  for (i = 0; i < pulseshape_size; i++)
    tlist[i] = t0+i*pulsewidth0/(pulseshape_size-1);
  acc = gsl_interp_accel_alloc();
  spline = gsl_spline_alloc(gsl_interp_cspline, pulseshape_size);
  gsl_spline_init(spline, tlist, pulseshape, pulseshape_size);
  free(tlist);
  //printf("test1\n");
  //Ramp Shape Interpolation
  tlist = (double *)malloc(sizeof(double)*rampshape_size);
  for (i = 0; i < rampshape_size; i++)
    tlist[i] = t0+rampon+i*(pulsewidth0-rampon-rampoff)/(rampshape_size-1);
  //printf("test2\n");
  acc1 = gsl_interp_accel_alloc();
  spline1 = gsl_spline_alloc(gsl_interp_cspline, rampshape_size);
  gsl_spline_init(spline1, tlist, rampshape, rampshape_size);
  //printf("test2\n");
  //Initialization
  psi = (struct complex *)malloc(sizeof(struct complex)*HSize*aoi_size);
  for (i = 0; i < HSize; i++)
  {
    psi[i].real = 0;
    psi[i].imag = 0;
  }
  psi[index].real = 1;
  Runge_Kutta(psi, t0, aoi, aoi_size, stepsize);
  for (i = 0; i < HSize*aoi_size; i++)
  {
    psi_real[i] = psi[i].real;
    psi_imag[i] = psi[i].imag;
  }
  free(psi);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  gsl_spline_free (spline1);
  gsl_interp_accel_free (acc1);
}
