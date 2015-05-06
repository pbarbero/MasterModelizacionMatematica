#define _GNU_SOURCE
#ifdef schedYes
#include <sched.h>
#endif
#include <stdio.h> 
#include <unistd.h> 
#include <omp.h> 

void affinity2_( ) {
#ifdef schedYes
  int NCPUs = sysconf(_SC_NPROCESSORS_CONF);
  cpu_set_t new_mask; 

  int tid = NCPUs/2 * omp_get_thread_num();
  CPU_ZERO( &new_mask );
  int i;
  for( i=0; i<NCPUs/2; i++ ) {
    CPU_SET( tid+i, &new_mask ); 
  }

  if (sched_setaffinity(0, sizeof(cpu_set_t), &new_mask) == -1) {
    printf("Error: sched_setaffinity(%d, sizeof(cpu_set_t), &new_mask)\n", tid); 
  } 
  /* printf("affinity2: tid=%d new_mask=%08X \n", tid, *(unsigned int*)(&new_mask)); */
#else
  if( !omp_get_thread_num() ) printf("Option --core-affinity cannot be used since the package is not currently compiled with --with-sched\n");
#endif 
}
