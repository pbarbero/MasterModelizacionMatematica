#define _GNU_SOURCE
#ifdef schedYes 
#include <sched.h>
#endif 
#include <stdio.h> 
#include <unistd.h> 
#include <omp.h> 

void affinity_( ) {
#ifdef schedYes 
  /* Total number of processors divided into two groups */
  int NCPUs = sysconf(_SC_NPROCESSORS_CONF) / 2;
  cpu_set_t new_mask; 
  cpu_set_t was_mask; 
  cpu_set_t mask; 

  int tid = omp_get_thread_num(); 
  CPU_ZERO(&new_mask);

  CPU_SET( tid%(2*NCPUs), &new_mask); 
  CPU_SET( (tid+NCPUs)%(2*NCPUs), &new_mask); 
  if (sched_getaffinity(0, sizeof(cpu_set_t), &was_mask) == -1) { 
    printf("Error: sched_getaffinity(%d, sizeof(cpu_set_t), &was_mask)\n", tid);
  } 
  CPU_AND( &mask, &new_mask, &was_mask );
  if (sched_setaffinity(0, sizeof(cpu_set_t), &mask) == -1) {
    printf("Error: sched_setaffinity(%d, sizeof(mask), &mask)\n", tid); 
  } 
  /* printf("affinity: tid=%d mask=%08X was_mask=%08X\n", tid, *(unsigned int*)(&mask), *(unsigned int*)(&was_mask)); */
#else
  /* if( !omp_get_thread_num() ) printf("It cannot used --core-affinity since the package is not currently compiled with --with-sched\n"); */
#endif 
}

