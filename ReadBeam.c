#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <math.h>

void get_beam_amp_phase(gsl_vector *beam_info, gsl_matrix *beam_amp, gsl_matrix *beam_phase, double x, double y, double *amp, double *phase)
{
  double beam_center_x = gsl_vector_get(beam_info,5);
  double beam_center_y = gsl_vector_get(beam_info,6);
  double dx = gsl_vector_get(beam_info,3);
  x = x/dx+beam_center_x;
  y = y/dx+beam_center_y;
  printf("%f\t%f\n",x,y);
  int j = floor(x);
  int i = floor(y);
  printf("%d\t%d\n",j,i);
  double amp0 = gsl_matrix_get(beam_amp,i,j)+(gsl_matrix_get(beam_amp,i,j+1)-gsl_matrix_get(beam_amp,i,j))*(x-j);
  double amp1 = gsl_matrix_get(beam_amp,i+1,j)+(gsl_matrix_get(beam_amp,i+1,j+1)-gsl_matrix_get(beam_amp,i+1,j))*(x-j);
  *amp = amp0+(amp1-amp0)*(y-i);
  i = i-1;
  j = j-1;
  double phi0 = gsl_matrix_get(beam_phase,i,j)+(gsl_matrix_get(beam_phase,i,j+1)-gsl_matrix_get(beam_phase,i,j))*(x-j);
  double phi1 = gsl_matrix_get(beam_phase,i+1,j)+(gsl_matrix_get(beam_phase,i+1,j+1)-gsl_matrix_get(beam_phase,i+1,j))*(x-j);
  *phase = phi0+(phi1-phi0)*(y-i);
}

gsl_matrix *beam_amp;
gsl_matrix *beam_phase;
gsl_vector *beam_info;

int
main (void)
{
  FILE *f_in_info = fopen("beam/beam_info.txt", "r");
  FILE *f_in_amp = fopen("beam/beam_amp.txt", "r");
  FILE *f_in_phase = fopen("beam/beam_phase.txt", "r");
  beam_info = gsl_vector_alloc(8);
  gsl_vector_fscanf(f_in_info, beam_info);
  int size_y = round(gsl_vector_get(beam_info,0));
  int size_x = round(gsl_vector_get(beam_info,1));
  beam_amp = gsl_matrix_alloc(size_y, size_x);
  beam_phase = gsl_matrix_alloc(size_y-2, size_x-2);
  gsl_matrix_fscanf(f_in_amp, beam_amp);
  gsl_matrix_fscanf(f_in_phase, beam_phase);
  double x, y;
  double *amp, *phase;
  amp = (double *)malloc(sizeof(double));
  phase = (double *)malloc(sizeof(double));
  scanf("%lf\t%lf", &x, &y);
  get_beam_amp_phase(beam_info, beam_amp, beam_phase, x, y, amp, phase);
  printf("amp=%.10lf\n",*amp);
  printf("phase=%.10lf\n",*phase);
  return 0;
}
