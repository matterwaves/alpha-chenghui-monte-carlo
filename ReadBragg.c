#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf.h>
#include <math.h>

gsl_interp_accel *acc0, *acc1;
gsl_spline *amp0, *amp1;

double contrast_u(double intensity)
{
  if (intensity < 0.1 || intensity > 1.25)
    return 0;
  double eff0 = gsl_spline_eval(amp0, intensity, acc0);
  double eff1 = gsl_spline_eval(amp1, intensity, acc1);
  eff0 = eff0*eff0;
  eff1 = eff1*eff1;
  return (4*eff1*eff1*eff0*eff0)/(eff1*eff1*eff1*eff0+eff1*eff0*eff0*eff0+2*eff1*eff1*eff0*eff0);
}

double contrast_l(double intensity)
{
  if (intensity < 0.1 || intensity > 1.25)
    return 0;
  double eff0 = gsl_spline_eval(amp0, intensity, acc0);
  double eff1 = gsl_spline_eval(amp1, intensity, acc1);
  eff0 = eff0*eff0;
  eff1 = eff1*eff1;
  return (4*eff1*eff1*eff0*eff0)/(eff1*eff1*eff1*eff1+eff0*eff0*eff0*eff0+eff1*eff1*eff1*eff0+eff1*eff0*eff0*eff0);
}

double Bloch_eff(double intensity)
{
  double *Bloch_Params = (double *)malloc(sizeof(double)*4);
  Bloch_Params[0] = 0.194;
  Bloch_Params[1] = 0.0179;
  Bloch_Params[2] = 0.215;
  Bloch_Params[3] = 322;
  intensity = intensity*1000;
  return Bloch_Params[0]*gsl_sf_erf(Bloch_Params[1]*(intensity-Bloch_Params[3]))+Bloch_Params[2];
}

int
main (void)
{
  int Bragg_size = 100;
  int i;
  gsl_matrix * Bragg_amp = gsl_matrix_alloc(Bragg_size, 3);
  double * intensity = (double *)malloc(sizeof(double)*Bragg_size);
  double * a0 = (double *)malloc(sizeof(double)*Bragg_size);
  double * a1 = (double *)malloc(sizeof(double)*Bragg_size);
  acc0 = gsl_interp_accel_alloc();
  amp0 = gsl_spline_alloc(gsl_interp_cspline, Bragg_size);
  acc1 = gsl_interp_accel_alloc();
  amp1 = gsl_spline_alloc(gsl_interp_cspline, Bragg_size);
  FILE *f_Bragg_in = fopen("Bragg_power.txt","r");
  gsl_matrix_fscanf(f_Bragg_in, Bragg_amp);
  for (i = 0; i < Bragg_size; i++)
  {
    intensity[i] = gsl_matrix_get(Bragg_amp, i, 0);
    a0[i] = sqrt(gsl_matrix_get(Bragg_amp, i, 1));
    a1[i] = sqrt(gsl_matrix_get(Bragg_amp, i, 2));
  }

  gsl_spline_init(amp0, intensity, a0, Bragg_size);
  gsl_spline_init(amp1, intensity, a1, Bragg_size);
  FILE *f_out = fopen("test.txt", "w");
  for (i = 0; i < Bragg_size; i++)
  {
    double intensity = i*0.012;
    fprintf(f_out, "%f\t%f\n", intensity, Bloch_eff(intensity));
  }
  return 0;
}
