#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_vector.h>
#include <time.h>
#include <math.h>


double Detuning, fL, c, k, hbar, sigma_x, sigma_v, kB, M, t0, T, Tp1, TB;
double Tp2, ztop, vr, omega_r, v_2nd_Bragg, g, Bragg_peak_rel, Bloch_peak_rel;
int n, N;

//Bragg Bloch efficiency
gsl_interp_accel *acc0, *acc1, *acc2;
gsl_spline *Bragg_amp0, *Bragg_amp1, *Bloch_eff;
int Bragg_size = 100, Bloch_size = 15;
gsl_matrix * Bragg_amp, * Bloch_data;

double contrast_u(double intensity)
{
  if (intensity < 0.1 || intensity > 1.25)
    return 0;
  double eff0 = gsl_spline_eval(Bragg_amp0, intensity, acc0);
  double eff1 = gsl_spline_eval(Bragg_amp1, intensity, acc1);
  eff0 = eff0*eff0;
  eff1 = eff1*eff1;
  return (4*eff1*eff1*eff0*eff0)/(eff1*eff1*eff1*eff0+eff1*eff0*eff0*eff0+2*eff1*eff1*eff0*eff0);
}

double contrast_l(double intensity)
{
  if (intensity < 0.1 || intensity > 1.25)
    return 0;
  double eff0 = gsl_spline_eval(Bragg_amp0, intensity, acc0);
  double eff1 = gsl_spline_eval(Bragg_amp1, intensity, acc1);
  eff0 = eff0*eff0;
  eff1 = eff1*eff1;
  return (4*eff1*eff1*eff0*eff0)/(eff1*eff1*eff1*eff1+eff0*eff0*eff0*eff0+eff1*eff1*eff1*eff0+eff1*eff0*eff0*eff0);
}

double Bloch_amp(double intensity)
{
  if (intensity < 0.1681038251290286)
    return 0;
  if (intensity > 0.22)
    return pow(0.728575225956,0.1);
  return pow(gsl_spline_eval(Bloch_eff, intensity, acc2),0.1);
}

void Read_Bragg_info()
{
  double * intensity = (double *)malloc(sizeof(double)*Bragg_size);
  double * a0 = (double *)malloc(sizeof(double)*Bragg_size);
  double * a1 = (double *)malloc(sizeof(double)*Bragg_size);
  double * intensity_B = (double *)malloc(sizeof(double)*Bloch_size);
  double * Bloch_eff0 = (double *)malloc(sizeof(double)*Bloch_size);
  int i;
  Bragg_amp = gsl_matrix_alloc(Bragg_size, 3);
  acc0 = gsl_interp_accel_alloc();
  Bragg_amp0 = gsl_spline_alloc(gsl_interp_cspline, Bragg_size);
  acc1 = gsl_interp_accel_alloc();
  Bragg_amp1 = gsl_spline_alloc(gsl_interp_cspline, Bragg_size);
  FILE *f_Bragg_in = fopen("Bragg_power.txt","r");
  gsl_matrix_fscanf(f_Bragg_in, Bragg_amp);
  for (i = 0; i < Bragg_size; i++)
  {
    intensity[i] = gsl_matrix_get(Bragg_amp, i, 0);
    a0[i] = sqrt(gsl_matrix_get(Bragg_amp, i, 1));
    a1[i] = sqrt(gsl_matrix_get(Bragg_amp, i, 2));
  }
  gsl_spline_init(Bragg_amp0, intensity, a0, Bragg_size);
  gsl_spline_init(Bragg_amp1, intensity, a1, Bragg_size);
  fclose(f_Bragg_in);

  Bloch_data = gsl_matrix_alloc(Bloch_size, 2);
  acc2 = gsl_interp_accel_alloc();
  Bloch_eff = gsl_spline_alloc(gsl_interp_cspline, Bloch_size);
  FILE *f_Bloch_in = fopen("Bloch_power.txt","r");
  gsl_matrix_fscanf(f_Bloch_in, Bloch_data);
  for (i = 0; i < Bloch_size; i++)
  {
    intensity_B[i] = gsl_matrix_get(Bloch_data, i, 0);
    Bloch_eff0[i] = gsl_matrix_get(Bloch_data, i, 1);
  }
  gsl_spline_init(Bloch_eff, intensity_B, Bloch_eff0, Bloch_size);
  fclose(f_Bragg_in);
}

//Beam amplitude and phase
gsl_matrix *beam_amp_up, *beam_amp_down;
gsl_matrix *beam_phase_up, *beam_phase_down;
gsl_vector *beam_info_up, *beam_info_down;

void get_beam_amp_phase(gsl_vector *beam_info, gsl_matrix *beam_amp, gsl_matrix *beam_phase, double x, double y, double *amp, double *phase)
{
  double beam_center_x = gsl_vector_get(beam_info,5);
  double beam_center_y = gsl_vector_get(beam_info,6);
  double dx = gsl_vector_get(beam_info,3);
  //printf("%f\t%f\n",x,y);
  x = x/dx+beam_center_x;
  y = y/dx+beam_center_y;
  //printf("%f\t%f\n",x,y);
  int j = floor(x);
  int i = floor(y);
  //printf("%d\t%d\n",j,i);
  //printf("test_end\n");
  int size_y = beam_amp->size1;
  int size_x = beam_amp->size2;
  if (i < 2 || j < 2 || i >= size_y-2 || j >= size_x-2)
  {
    *amp = 0;
    *phase = 0;
  }
  else
  {
    double amp0 = gsl_matrix_get(beam_amp,i,j)+(gsl_matrix_get(beam_amp,i,j+1)-gsl_matrix_get(beam_amp,i,j))*(x-j);
    double amp1 = gsl_matrix_get(beam_amp,i+1,j)+(gsl_matrix_get(beam_amp,i+1,j+1)-gsl_matrix_get(beam_amp,i+1,j))*(x-j);
    *amp = amp0+(amp1-amp0)*(y-i);
    i = i-1;
    j = j-1;
    x = x-1;
    y = y-1;
    double phi0 = gsl_matrix_get(beam_phase,i,j)+(gsl_matrix_get(beam_phase,i,j+1)-gsl_matrix_get(beam_phase,i,j))*(x-j);
    double phi1 = gsl_matrix_get(beam_phase,i+1,j)+(gsl_matrix_get(beam_phase,i+1,j+1)-gsl_matrix_get(beam_phase,i+1,j))*(x-j);
    *phase = phi0+(phi1-phi0)*(y-i);
  }
  //printf("amp, phase = %.15f, %.15f\n", *amp, *phase);
}

void Read_beam_data()
{
  //Beam Up
  FILE *f_in_info = fopen("beam_up/beam_info.txt", "r");
  FILE *f_in_amp = fopen("beam_up/beam_amp.txt", "r");
  FILE *f_in_phase = fopen("beam_up/beam_phase.txt", "r");
  beam_info_up = gsl_vector_alloc(8);
  beam_info_down = gsl_vector_alloc(8);
  gsl_vector_fscanf(f_in_info, beam_info_up);
  int size_y = round(gsl_vector_get(beam_info_up,0));
  int size_x = round(gsl_vector_get(beam_info_up,1));
  beam_amp_up = gsl_matrix_alloc(size_y, size_x);
  beam_phase_up = gsl_matrix_alloc(size_y-2, size_x-2);
  gsl_matrix_fscanf(f_in_amp, beam_amp_up);
  gsl_matrix_fscanf(f_in_phase, beam_phase_up);
  fclose(f_in_info);
  fclose(f_in_amp);
  fclose(f_in_phase);
  //Beam Down
  f_in_info = fopen("beam_down/beam_info.txt", "r");
  f_in_amp = fopen("beam_down/beam_amp.txt", "r");
  f_in_phase = fopen("beam_down/beam_phase.txt", "r");
  gsl_vector_fscanf(f_in_info, beam_info_down);
  size_y = round(gsl_vector_get(beam_info_down,0));
  size_x = round(gsl_vector_get(beam_info_down,1));
  beam_amp_down = gsl_matrix_alloc(size_y, size_x);
  beam_phase_down = gsl_matrix_alloc(size_y-2, size_x-2);
  gsl_matrix_fscanf(f_in_amp, beam_amp_down);
  gsl_matrix_fscanf(f_in_phase, beam_phase_down);
}


//Parameters
//Laser
void set_parameters(void)
{
  Detuning = 14E+9;
  fL =  3.517309021e+14+90E+6+80E+6+Detuning;
  c = 2.99792458e+8;
  k = 2*M_PI*fL/c;
  hbar = 1.054571800E-34;
  //Atom distribution
  sigma_x = 0.002;
  sigma_v = 0.0035;
  kB = 1.38064852e-23;
  M = 2.20694650e-25;
  //Atom Interferometer Parameters
  t0 = 1.24-1.12;
  T = 0.060;
  Tp1 = 0.005;
  TB = 0.008;
  Tp2 = 0.015;
  ztop = 65*0.0254;
  vr = hbar*k/M;
  omega_r = hbar*k*k/2/M;
  n = 5;
  N = 125;
  //Atom trajectory
  v_2nd_Bragg = 2.0480821101976;
  g = 9.79958;
  //Parameters
  Bragg_peak_rel = 1.05;
  Bloch_peak_rel = 0.2;
}


void move(double *x, double *y, gsl_vector *z, double *vx, double *vy, gsl_vector *vz, double t)
{
  *x = *x + *vx * t;
  *y = *y + *vy * t;
  gsl_vector_set(z, 0, gsl_vector_get(z, 0)+gsl_vector_get(vz, 0)*t-0.5*g*t*t);
  gsl_vector_set(z, 1, gsl_vector_get(z, 1)+gsl_vector_get(vz, 1)*t-0.5*g*t*t);
  gsl_vector_set(vz, 0, gsl_vector_get(vz, 0)-g*t);
  gsl_vector_set(vz, 1, gsl_vector_get(vz, 1)-g*t);
}

double mean(gsl_vector *data)
{
  int size = data->size;
  double sum = 0;
  int i;
  for (i = 0; i < size; i++)
  {
    sum += gsl_vector_get(data, i);
  }
  return sum/size;
}

gsl_vector *xset, *yset, *vxset, *vyset, *Bloch_eff_stat;
gsl_matrix *phase, *amp;
gsl_matrix *peak;
const gsl_rng_type * Type;
gsl_rng * r;
FILE *f_out;

int MonteCarlo(int atom_number, int bin_size, int seed)
{
  int i, max_index;
  double x, y, vx, vy;
  double amp_up, amp_down, intensity, phase_up, phase_down;
  double max_amp = 0;
  gsl_vector *z = gsl_vector_calloc(2);
  gsl_vector *vz = gsl_vector_alloc(2);
  gsl_vector *current_amp;
  clock_t begin = clock();
  xset = gsl_vector_alloc(atom_number);
  yset = gsl_vector_alloc(atom_number);
  vxset = gsl_vector_alloc(atom_number);
  vyset = gsl_vector_alloc(atom_number);
  Bloch_eff_stat = gsl_vector_calloc(atom_number);
  gsl_vector_set_all(Bloch_eff_stat, 1);
  phase = gsl_matrix_calloc(atom_number, 2);
  amp = gsl_matrix_alloc(atom_number, 4);
  current_amp = gsl_vector_calloc(4);
  peak = gsl_matrix_calloc(bin_size*2, 4);
  gsl_matrix_set_all(amp, 1);
  double vz0 = v_2nd_Bragg+g*(1.32-1.12);
  gsl_vector_set_all(vz, vz0);
  /* create a generator chosen by the
  environment variable GSL_RNG_TYPE */
  gsl_rng_env_setup();
  Type = gsl_rng_default;
  r = gsl_rng_alloc (Type);
  gsl_rng_set(r, seed);
  /* print n random variates chosen from
  the poisson distribution with mean
  parameter mu */
  //printf("generator type: %s\n", gsl_rng_name(r));
  //printf ("seed = %lu\n", gsl_rng_default_seed);
  for (i = 0; i < atom_number; i++)
  {
    //if (i%10000 == 0)
      //printf("Atom#%d\n",i);
    gsl_ran_bivariate_gaussian (r, sigma_x, sigma_x, 0, &x, &y);
    gsl_ran_bivariate_gaussian (r, sigma_v, sigma_v, 0, &vx, &vy);
    //vx = 0;
    //vy = 0;
    //x = 0.001;
    //y = 0;
    //printf("x = %.10f, y = %.10f, z = %.10f, vx = %.10f, vy = %.10f, vz = %.10f\n", x, y, gsl_vector_get(z,0), vx, vy, gsl_vector_get(vz,0));
    move(&x, &y, z, &vx, &vy, vz, 1.32-T-1.12);
    //First pulse
    //t = -T
    //printf("amp_i, phase_i = %.10f, %.10f\n", gsl_matrix_get(amp, i, 0), gsl_matrix_get(phase, i, 0));
    //printf("x = %.10f, y = %.10f, z = %.10f, vx = %.10f, vy = %.10f, vz = %.10f\n", x, y, gsl_vector_get(z,0), vx, vy, gsl_vector_get(vz,0));
    get_beam_amp_phase(beam_info_up, beam_amp_up, beam_phase_up, x, y, &amp_up, &phase_up);
    get_beam_amp_phase(beam_info_down, beam_amp_down, beam_phase_down, x, y, &amp_down, &phase_down);
    intensity = Bragg_peak_rel*amp_up*amp_down;
    //printf("amp_up, phase_up = %.15f, %.15f\n", amp_up, phase_up);
    //#print(amp_up, amp_down)
    gsl_matrix_set(amp, i, 0, gsl_spline_eval(Bragg_amp1, intensity, acc1));
    gsl_matrix_set(amp, i, 1, gsl_spline_eval(Bragg_amp0, intensity, acc0));
    //printf("intensity = %.10f\n", intensity);
    //printf("amp0, amp1 = %.10f, %.10f\n", gsl_matrix_get(amp, i, 0), gsl_matrix_get(amp, i, 1));
    //printf("Bragg_amp0 = %.10f\n", gsl_spline_eval(Bragg_amp0, intensity, acc0));
    //#print(amp[i,:])
    move(&x, &y, z, &vx, &vy, vz, T);
    //printf("amp_i, phase_i = %.10f, %.10f\n", gsl_matrix_get(amp, i, 0), gsl_matrix_get(phase, i, 0));
    //Second pulse
    //t = 0
    //printf("x = %.10f, y = %.10f, z = %.10f, vx = %.10f, vy = %.10f, vz = %.10f\n", x, y, gsl_vector_get(z,0), vx, vy, gsl_vector_get(vz,0));
    get_beam_amp_phase(beam_info_up, beam_amp_up, beam_phase_up, x, y, &amp_up, &phase_up);
    get_beam_amp_phase(beam_info_down, beam_amp_down, beam_phase_down, x, y, &amp_down, &phase_down);
    intensity = Bragg_peak_rel*amp_up*amp_down;
    gsl_matrix_get_row(current_amp, amp, i);
    gsl_matrix_set(amp, i, 0, gsl_vector_get(current_amp, 0)*gsl_spline_eval(Bragg_amp1, intensity, acc1));
    gsl_matrix_set(amp, i, 1, gsl_vector_get(current_amp, 1)*gsl_spline_eval(Bragg_amp0, intensity, acc0));
    gsl_matrix_set(amp, i, 2, gsl_vector_get(current_amp, 0)*gsl_spline_eval(Bragg_amp0, intensity, acc0));
    gsl_matrix_set(amp, i, 3, gsl_vector_get(current_amp, 1)*gsl_spline_eval(Bragg_amp1, intensity, acc1));
    //printf("intensity = %.10f\n", intensity);
    //printf("amp0, amp1, amp2, amp3 = %.10f, %.10f, %.10f, %.10f\n", gsl_matrix_get(amp, i, 0), gsl_matrix_get(amp, i, 1), gsl_matrix_get(amp, i, 2), gsl_matrix_get(amp, i, 3));
    //printf("Bragg_amp0 = %.10f\n", gsl_spline_eval(Bragg_amp0, intensity, acc0));
    //printf("x = %f, y = %f, z = %f, vx = %f, vy = %f, vz = %f\n", x, y, gsl_vector_get(z,0), vx, vy, gsl_vector_get(vz,0));
    gsl_matrix_set(phase, i, 0, gsl_matrix_get(phase, i, 0)+2*n*n*omega_r*T*(phase_up+phase_down));
    gsl_matrix_set(phase, i, 1, gsl_matrix_get(phase, i, 1)-2*n*n*omega_r*T*(phase_up+phase_down));
    move(&x, &y, z, &vx, &vy, vz, Tp1);
    //printf("amp_i, phase_i = %.10f, %.10f\n", gsl_matrix_get(amp, i, 0), gsl_matrix_get(phase, i, 0));
    //#Bloch oscillations
    //t = Tp1
    int j;
    for (j = 0; j < 5; j++)
    {
      //printf("x = %.10f, y = %.10f, z = %.10f, vx = %.10f, vy = %.10f, vz = %.10f\n", x, y, gsl_vector_get(z,0), vx, vy, gsl_vector_get(vz,0));
      get_beam_amp_phase(beam_info_up, beam_amp_up, beam_phase_up, x, y, &amp_up, &phase_up);
      get_beam_amp_phase(beam_info_down, beam_amp_down, beam_phase_down, x, y, &amp_down, &phase_down);
      intensity = Bloch_peak_rel*amp_up*amp_down;
      gsl_vector_set(Bloch_eff_stat, i, gsl_vector_get(Bloch_eff_stat, i)*Bloch_amp(intensity)*Bloch_amp(intensity));
      gsl_matrix_set(amp, i, 0, gsl_matrix_get(amp, i, 0)*Bloch_amp(intensity));
      gsl_matrix_set(amp, i, 1, gsl_matrix_get(amp, i, 1)*Bloch_amp(intensity));
      gsl_matrix_set(amp, i, 2, gsl_matrix_get(amp, i, 2)*Bloch_amp(intensity));
      gsl_matrix_set(amp, i, 3, gsl_matrix_get(amp, i, 3)*Bloch_amp(intensity));
      //printf("intensity = %.10f\n", intensity);
      //#print(amp[i,:])
      gsl_matrix_set(phase, i, 0, gsl_matrix_get(phase, i, 0)+4*n*(N/5)*omega_r*T*(phase_up+phase_down));
      gsl_matrix_set(phase, i, 1, gsl_matrix_get(phase, i, 1)-4*n*(N/5)*omega_r*T*(phase_up+phase_down));
      gsl_vector_set(vz, 0, gsl_vector_get(vz, 0)-2*(N/5)*vr);
      gsl_vector_set(vz, 1, gsl_vector_get(vz, 1)+2*(N/5)*vr);
      move(&x, &y, z, &vx, &vy, vz, TB/5);
    }
    move(&x, &y, z, &vx, &vy, vz, Tp2-TB);
    //printf("amp_i, phase_i = %.10f, %.10f\n", gsl_matrix_get(amp, i, 0), gsl_matrix_get(phase, i, 0));
    //Third pulse
    //t = Tp1+Tp2
    //printf("x = %.10f, y = %.10f, z = %.10f, vx = %.10f, vy = %.10f, vz = %.10f\n", x, y, gsl_vector_get(z,0), vx, vy, gsl_vector_get(vz,0));
    get_beam_amp_phase(beam_info_up, beam_amp_up, beam_phase_up, x, y, &amp_up, &phase_up);
    get_beam_amp_phase(beam_info_down, beam_amp_down, beam_phase_down, x, y, &amp_down, &phase_down);
    intensity = Bragg_peak_rel*amp_up*amp_down;
    //printf("Bragg_amp0 = %.10f\n", gsl_spline_eval(Bragg_amp0, intensity, acc0));
    gsl_matrix_set(amp, i, 0, gsl_matrix_get(amp, i, 0)*gsl_spline_eval(Bragg_amp1, intensity, acc1));
    gsl_matrix_set(amp, i, 1, gsl_matrix_get(amp, i, 1)*gsl_spline_eval(Bragg_amp0, intensity, acc0));
    gsl_matrix_set(amp, i, 2, gsl_matrix_get(amp, i, 2)*gsl_spline_eval(Bragg_amp0, intensity, acc0));
    gsl_matrix_set(amp, i, 3, gsl_matrix_get(amp, i, 3)*gsl_spline_eval(Bragg_amp1, intensity, acc1));
    gsl_matrix_set(phase, i, 0, gsl_matrix_get(phase, i, 0)+2*n*n*omega_r*T*(phase_up+phase_down));
    gsl_matrix_set(phase, i, 1, gsl_matrix_get(phase, i, 1)-2*n*n*omega_r*T*(phase_up+phase_down));
    move(&x, &y, z, &vx, &vy, vz, T);
    //printf("intensity = %.10f\n", intensity);
    //printf("amp0, amp1, amp2, amp3 = %.10f, %.10f, %.10f, %.10f\n", gsl_matrix_get(amp, i, 0), gsl_matrix_get(amp, i, 1), gsl_matrix_get(amp, i, 2), gsl_matrix_get(amp, i, 3));
    //printf("amp_i, phase_i = %.10f, %.10f\n", gsl_matrix_get(amp, i, 0), gsl_matrix_get(phase, i, 0));
    //#Fourth pulse
    //t = Tp1+Tp2+T
    //printf("x = %.10f, y = %.10f, z = %.10f, vx = %.10f, vy = %.10f, vz = %.10f\n", x, y, gsl_vector_get(z,0), vx, vy, gsl_vector_get(vz,0));
    get_beam_amp_phase(beam_info_up, beam_amp_up, beam_phase_up, x, y, &amp_up, &phase_up);
    get_beam_amp_phase(beam_info_down, beam_amp_down, beam_phase_down, x, y, &amp_down, &phase_down);
    intensity = Bragg_peak_rel*amp_up*amp_down;
    //printf("intensity = %.10f\n", intensity);
    //printf("Bragg_amp0 = %.10f\n", gsl_spline_eval(Bragg_amp0, intensity, acc0));
    gsl_matrix_set(phase, i, 0, gsl_matrix_get(phase, i, 0)+(4*n*n+4*n*N)*omega_r*T*(phase_up+phase_down));
    gsl_matrix_set(phase, i, 1, gsl_matrix_get(phase, i, 1)-(4*n*n+4*n*N)*omega_r*T*(phase_up+phase_down));
    for (j = 0; j < 30; j++)
    {
      gsl_vector *current_signal;
      current_signal = gsl_vector_calloc(4);
      double common_phase, phaseu, phasel;
      phasel = gsl_matrix_get(phase, i, 0)-M_PI/2;
      phaseu = gsl_matrix_get(phase, i, 1);
      //phasel = -M_PI/2;
      //phaseu = 0.001;
      common_phase = j/30.0*2*M_PI;
      //printf("common_phase = %.10e\n",common_phase);
      //printf("phasel = %.10e\n",phasel);
      //printf("phaseu = %.10e\n",phaseu);
      gsl_vector_set(current_signal, 0, gsl_matrix_get(amp, i, 0)*gsl_matrix_get(amp, i, 0)*gsl_spline_eval(Bragg_amp0, intensity, acc0)*gsl_spline_eval(Bragg_amp0, intensity, acc0)
      +gsl_matrix_get(amp, i, 1)*gsl_matrix_get(amp, i, 1)*gsl_spline_eval(Bragg_amp1, intensity, acc1)*gsl_spline_eval(Bragg_amp1, intensity, acc1)
      +2*gsl_matrix_get(amp, i, 0)*gsl_matrix_get(amp, i, 1)*gsl_spline_eval(Bragg_amp0, intensity, acc0)*gsl_spline_eval(Bragg_amp1, intensity, acc1)*gsl_sf_sin(phasel+common_phase));

      gsl_vector_set(current_signal, 1, gsl_matrix_get(amp, i, 0)*gsl_matrix_get(amp, i, 0)*gsl_spline_eval(Bragg_amp1, intensity, acc1)*gsl_spline_eval(Bragg_amp1, intensity, acc1)
      +gsl_matrix_get(amp, i, 1)*gsl_matrix_get(amp, i, 1)*gsl_spline_eval(Bragg_amp0, intensity, acc0)*gsl_spline_eval(Bragg_amp0, intensity, acc0)
      -2*gsl_matrix_get(amp, i, 0)*gsl_matrix_get(amp, i, 1)*gsl_spline_eval(Bragg_amp0, intensity, acc0)*gsl_spline_eval(Bragg_amp1, intensity, acc1)*gsl_sf_sin(phasel+common_phase));

      gsl_vector_set(current_signal, 2, gsl_matrix_get(amp, i, 2)*gsl_matrix_get(amp, i, 2)*gsl_spline_eval(Bragg_amp0, intensity, acc0)*gsl_spline_eval(Bragg_amp0, intensity, acc0)
      +gsl_matrix_get(amp, i, 3)*gsl_matrix_get(amp, i, 3)*gsl_spline_eval(Bragg_amp1, intensity, acc1)*gsl_spline_eval(Bragg_amp1, intensity, acc1)
      +2*gsl_matrix_get(amp, i, 2)*gsl_matrix_get(amp, i, 3)*gsl_spline_eval(Bragg_amp0, intensity, acc0)*gsl_spline_eval(Bragg_amp1, intensity, acc1)*gsl_sf_sin(phaseu+common_phase));

      gsl_vector_set(current_signal, 3, gsl_matrix_get(amp, i, 2)*gsl_matrix_get(amp, i, 2)*gsl_spline_eval(Bragg_amp1, intensity, acc1)*gsl_spline_eval(Bragg_amp1, intensity, acc1)
      +gsl_matrix_get(amp, i, 3)*gsl_matrix_get(amp, i, 3)*gsl_spline_eval(Bragg_amp0, intensity, acc0)*gsl_spline_eval(Bragg_amp0, intensity, acc0)
      -2*gsl_matrix_get(amp, i, 2)*gsl_matrix_get(amp, i, 3)*gsl_spline_eval(Bragg_amp0, intensity, acc0)*gsl_spline_eval(Bragg_amp1, intensity, acc1)*gsl_sf_sin(phaseu+common_phase));

      gsl_matrix_set(peak, j, 0, gsl_matrix_get(peak, j, 0)+gsl_vector_get(current_signal,0));
      gsl_matrix_set(peak, j, 1, gsl_matrix_get(peak, j, 1)+gsl_vector_get(current_signal,1));
      gsl_matrix_set(peak, j, 2, gsl_matrix_get(peak, j, 2)+gsl_vector_get(current_signal,2));
      gsl_matrix_set(peak, j, 3, gsl_matrix_get(peak, j, 3)+gsl_vector_get(current_signal,3));

      //printf("%.10e\t%.10e\n", (gsl_vector_get(current_signal,0)-gsl_vector_get(current_signal,1))/(gsl_vector_get(current_signal,0)+gsl_vector_get(current_signal,1)), (gsl_vector_get(current_signal,2)-gsl_vector_get(current_signal,3))/(gsl_vector_get(current_signal,2)+gsl_vector_get(current_signal,3)));
      //printf("%d = \n", j);
      //printf("%.10e, %.10e, %.10e, %.10e\n", gsl_matrix_get(peak, j, 0), gsl_matrix_get(peak, j, 1), gsl_matrix_get(peak, j, 2), gsl_matrix_get(peak, j, 3));
      gsl_vector_free(current_signal);
    }

    for (j = 0; j < 30; j++)
    {
      gsl_vector *current_signal;
      current_signal = gsl_vector_calloc(4);
      double common_phase, phaseu, phasel;
      phasel = gsl_matrix_get(phase, i, 0)+M_PI/2;
      phaseu = gsl_matrix_get(phase, i, 1);
      //phasel = M_PI/2;
      //phaseu = 0.001;
      common_phase = j/30.0*2*M_PI;
      gsl_vector_set(current_signal, 0, gsl_matrix_get(amp, i, 0)*gsl_matrix_get(amp, i, 0)*gsl_spline_eval(Bragg_amp0, intensity, acc0)*gsl_spline_eval(Bragg_amp0, intensity, acc0)
      +gsl_matrix_get(amp, i, 1)*gsl_matrix_get(amp, i, 1)*gsl_spline_eval(Bragg_amp1, intensity, acc1)*gsl_spline_eval(Bragg_amp1, intensity, acc1)
      +2*gsl_matrix_get(amp, i, 0)*gsl_matrix_get(amp, i, 1)*gsl_spline_eval(Bragg_amp0, intensity, acc0)*gsl_spline_eval(Bragg_amp1, intensity, acc1)*gsl_sf_sin(phasel+common_phase));

      gsl_vector_set(current_signal, 1, gsl_matrix_get(amp, i, 0)*gsl_matrix_get(amp, i, 0)*gsl_spline_eval(Bragg_amp1, intensity, acc1)*gsl_spline_eval(Bragg_amp1, intensity, acc1)
      +gsl_matrix_get(amp, i, 1)*gsl_matrix_get(amp, i, 1)*gsl_spline_eval(Bragg_amp0, intensity, acc0)*gsl_spline_eval(Bragg_amp0, intensity, acc0)
      -2*gsl_matrix_get(amp, i, 0)*gsl_matrix_get(amp, i, 1)*gsl_spline_eval(Bragg_amp0, intensity, acc0)*gsl_spline_eval(Bragg_amp1, intensity, acc1)*gsl_sf_sin(phasel+common_phase));

      gsl_vector_set(current_signal, 2, gsl_matrix_get(amp, i, 2)*gsl_matrix_get(amp, i, 2)*gsl_spline_eval(Bragg_amp0, intensity, acc0)*gsl_spline_eval(Bragg_amp0, intensity, acc0)
      +gsl_matrix_get(amp, i, 3)*gsl_matrix_get(amp, i, 3)*gsl_spline_eval(Bragg_amp1, intensity, acc1)*gsl_spline_eval(Bragg_amp1, intensity, acc1)
      +2*gsl_matrix_get(amp, i, 2)*gsl_matrix_get(amp, i, 3)*gsl_spline_eval(Bragg_amp0, intensity, acc0)*gsl_spline_eval(Bragg_amp1, intensity, acc1)*gsl_sf_sin(phaseu+common_phase));

      gsl_vector_set(current_signal, 3, gsl_matrix_get(amp, i, 2)*gsl_matrix_get(amp, i, 2)*gsl_spline_eval(Bragg_amp1, intensity, acc1)*gsl_spline_eval(Bragg_amp1, intensity, acc1)
      +gsl_matrix_get(amp, i, 3)*gsl_matrix_get(amp, i, 3)*gsl_spline_eval(Bragg_amp0, intensity, acc0)*gsl_spline_eval(Bragg_amp0, intensity, acc0)
      -2*gsl_matrix_get(amp, i, 2)*gsl_matrix_get(amp, i, 3)*gsl_spline_eval(Bragg_amp0, intensity, acc0)*gsl_spline_eval(Bragg_amp1, intensity, acc1)*gsl_sf_sin(phaseu+common_phase));

      gsl_matrix_set(peak, j+bin_size, 0, gsl_matrix_get(peak, j+bin_size, 0)+gsl_vector_get(current_signal,0));
      gsl_matrix_set(peak, j+bin_size, 1, gsl_matrix_get(peak, j+bin_size, 1)+gsl_vector_get(current_signal,1));
      gsl_matrix_set(peak, j+bin_size, 2, gsl_matrix_get(peak, j+bin_size, 2)+gsl_vector_get(current_signal,2));
      gsl_matrix_set(peak, j+bin_size, 3, gsl_matrix_get(peak, j+bin_size, 3)+gsl_vector_get(current_signal,3));
      //printf("%d = \n", j+bin_size);
      //printf("%.10e, %.10e, %.10e, %.10e\n", gsl_matrix_get(peak, j, 0), gsl_matrix_get(peak, j, 1), gsl_matrix_get(peak, j, 2), gsl_matrix_get(peak, j, 3));
      gsl_vector_free(current_signal);
    }

    move(&x, &y, z, &vx, &vy, vz, 1.93-1.32-T-Tp1-Tp2);
    gsl_vector_set(xset, i, x);
    gsl_vector_set(yset, i, y);
    gsl_vector_set(vxset, i, vx);
    gsl_vector_set(vyset, i, vy);
    //printf("amp0 = %.10f , phase0 = %.10f \n", gsl_matrix_get(amp, i, 0), gsl_matrix_get(phase, i, 0));
  }
  //printf("max_amp, max_index, phase, x, y, vx, vy = %.10lf, %d, %.10lf, %.10lf, %.10lf, %.10lf, %.10lf\n", max_amp, max_index, gsl_matrix_get(phase, max_index, 0), gsl_vector_get(xset, max_index), gsl_vector_get(yset, max_index), gsl_vector_get(vxset, max_index), gsl_vector_get(vyset, max_index));
  gsl_matrix_fprintf(f_out, peak, "%.10e");
  //gsl_vector *phase_sample, *amp_sample;
  //phase_sample = gsl_vector_alloc(atom_number);
  //amp_sample = gsl_vector_alloc(atom_number);
  //gsl_matrix_get_col(amp_sample, amp, 0);
  //gsl_matrix_get_col(phase_sample, phase, 0);
  //double mean_phase0 = mean(phase_sample, amp_sample)/16/n/(n+N);
  //gsl_matrix_get_col(amp_sample, amp, 1);
  //gsl_matrix_get_col(phase_sample, phase, 1);
  //double mean_phase1 = mean(phase_sample, amp_sample)/16/n/(n+N);
  //printf("Current phase = %.10f\n", mean_phase0+mean_phase1);
  //FILE *f_phase_out, *f_amp_out;
  //f_phase_out = fopen("phase_out.txt", "w");
  //f_amp_out = fopen("amp_out.txt", "w");
  //gsl_vector_fprintf(f_phase_out, phase_sample, "%.10lf");
  //gsl_vector_fprintf(f_amp_out, amp_sample, "%.10lf");
  //fclose(f_phase_out);
  //fclose(f_amp_out);
  gsl_rng_free (r);
  //gsl_vector_free(phase_sample);
  //gsl_vector_free(amp_sample);
  gsl_matrix_free(amp);
  gsl_vector_free(current_amp);
  gsl_matrix_free(peak);
  gsl_matrix_free(phase);
  gsl_vector_free(xset);
  gsl_vector_free(yset);
  gsl_vector_free(vxset);
  gsl_vector_free(vyset);
  clock_t end = clock();
  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  printf("Time spent = %f\n", time_spent);
  return 0;
}

int main()
{
  int i, repeats = 200, atom_number = 500000, bin_size = 30;
  double x, y, amp, phase;
  set_parameters();
  Read_Bragg_info();
  Read_beam_data();
  f_out = fopen("Results.txt", "w");
  fclose(f_out);
  f_out = fopen("Results.txt", "a");
  for (i = 0; i < repeats; i++)
  {
    Bragg_peak_rel = 1.2;
    Bloch_peak_rel = 0.24;
    printf("Atom#%d\n", i+1);
    int time_now = (int)time(NULL);
    printf("seed = %d\n",time_now+i);
    MonteCarlo(atom_number, bin_size, time_now+i);
    printf("Mean Bloch eff = %.10lf\n", mean(Bloch_eff_stat));
  }
  fclose(f_out);
  gsl_spline_free (Bragg_amp0);
  gsl_interp_accel_free (acc0);
  gsl_spline_free (Bragg_amp1);
  gsl_interp_accel_free (acc1);
  gsl_spline_free (Bloch_eff);
  gsl_interp_accel_free (acc2);
  gsl_matrix_free(Bragg_amp);
  gsl_matrix_free(Bloch_data);
  gsl_matrix_free(beam_amp_up);
  gsl_matrix_free(beam_amp_down);
  gsl_matrix_free(beam_phase_up);
  gsl_matrix_free(beam_phase_down);
  gsl_vector_free(beam_info_up);
  gsl_vector_free(beam_info_down);
}
