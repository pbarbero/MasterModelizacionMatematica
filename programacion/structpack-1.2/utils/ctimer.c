/*
 * Input:  none
 * Output: elapsed, milisegundos de tiempo real transcurridos desde la ultima vez que
 *                  fue llamada con gtime = 0.0
 *         ucpu,    milisegundos consumidos de cpu desde el inicio del programa
 *         scpu,    milisegundos consumidos de cpu por el sistema desde el inicio del programa
 *         gtime,   si 0.0 inicializa el contador elapsed
 * 
 * La rutina debe ejecutarse al menos dos veces. La primera con gtime = 0.0 pone en
 * marcha el cronómetro y devuelve 0 en elapsed. Los parámetros ucpu y scpu 
 * devuelven el valor correcto. La segunda, con gtime distinto de 0.0 devuelve 
 * en elapsed el tiempo real transcurrido entre la primera llamada y ésta.
 * La tercera y siguientes llamadas con gtime distinto de 0.0 devuelven el tiempo
 * transcurrido entre la ultima con gtime=0.0 y esta.
 *
 * Example:
 *
 *  double elapsed, ucpu, scpu, gtime=0.0;
 *
 *  ctimer(&elapsed, &ucpu, &scpu, &gtime);
 *  for (i=0;i<n;i++) {
 *     for (j=0;j<n;j++) {
 *       x[i] = A[i][j] * b[j];
 *     }
 *  }
 *  ctimer(&elapsed, &ucpu, &scpu, &gtime);
 *  printf("El producto x=A*b ha necesitado %f milisegundos\n", elapsed);
 *
*/
#include <time.h>
#include <sys/time.h>
#include <sys/times.h>
#include <unistd.h>
#include <cutils.h>

void ctimer_( double *elapsed, double *ucpu, double *scpu, double *gtime ) {
  struct timeval tm;
  struct tms sistema;

  struct timezone_obsoleto 
  {
     int   tz_minuteswest;
     int   tz_dsttime;
  } tz;

  double usegs;

  gettimeofday(&tm, &tz);
  times(&sistema);

  usegs = tm.tv_usec+tm.tv_sec*1E6;

  if (*gtime > 0.0)  {
    *elapsed = usegs - *gtime;
  } else {
    *elapsed = 0.0;
    *gtime = usegs;
  }

  *elapsed = *elapsed/1E6;
  *ucpu = (double)sistema.tms_utime/(double)CLOCKS_PER_SEC*1E4;
  *scpu = (double)sistema.tms_stime/(double)CLOCKS_PER_SEC*1E4;

} /* end ctimer */
