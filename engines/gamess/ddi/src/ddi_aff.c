/* -------------------------------------------------------------------- *\
 * Distributed Data Interface
 * ==========================
 *
 * Set CPU affinity for sockets
 * (MPI is assumed to have separate routines in MPI itself);
 * so that no MPI code is provided here to set the affinity.
 * The subroutines used below (sched_setaffinity etc) are inherent
 * to Linux with cores after 2.5.8; they may not be available for
 * other non-Linux UNIX...
 *
 * Should there be any linking problem, just set SOC_AFF to "no" in compddi.
 *
\* -------------------------------------------------------------------- */

 #define _GNU_SOURCE
 #include "ddi_base.h"

 # if defined DDI_SOC && defined SOC_AFF

 #include <sched.h>
 #include <unistd.h>
 #include <errno.h>

 void socket_aff(int np, int me, int my) {
 cpu_set_t set;
 int cpuid,npp,mee,err;

 /* fprintf(stderr," rank %d node %d.\n",me,my); */

 if(me<np)
   /* compute process affinity */
   cpuid=gv(ddinodes)[my].affinity[0];
 else
   /* data server affinity */
   cpuid=gv(ddinodes)[my].affinity[1];

/* 
   There are some special numbers for the affinity:
   -1 do not set the affinity 
   -2 set the affinity to be the core ID

    Usage: ** pay attention **
    core affinity is set in arguments to ddikick, 
    e.g., node1:cpus=3:aff=-2 or node1:cpus=1:aff=2  
    (a) For multi-core runs (e.g., node1:cpus=3), -2 is the only way.
        It means that a core affinity is set internally for every core
        according to its rank.
    (b) Logical nodes such as node1:cpus=1 node1:cpus=1 should set individual 
        values (do NOT use -2 for logical nodes). 
        It means that in "node1:cpus=1:aff=2" you request that this node1
        is run on CPU core #2. You can repeat as "node1:cpus=1:aff=3",
        asking for another node1 with 1 core to be run on core 3.
    Nota bene: do not use the same physical node with multiple logical nodes.
    Examples:
    1. allowed:
       node1:cpus=4:aff=-2 node2:cpus=4:aff=-2
       Why allowed: requested automatic core affinity on different nodes.
       node1:cpus=1:aff=0 node1:cpus=1:aff=1 node2:cpus=1:aff=0
       Why allowed: requested manual core affinity on the same nodes.
    2. Not allowed:
       node1:cpus=4:aff=2
       Why not allowed: manual core affinity (2) on a multicore node.
       node1:cpus=2:aff=-2 node1:cpus=2:aff=-2
       Why not allowed: automatic core affinity (-2) on a redundant node.

    The default is not to set affinity.
*/

 if(cpuid==-2)
 {
   DDI_SMP_NProc(&npp,&mee);
   if(me<np)
     /* compute process affinity */
     cpuid=mee;
   else
     /* data server affinity */
     cpuid=mee+npp;
 }

 if(cpuid>=0)
 {
   CPU_ZERO(&set);
   CPU_SET(cpuid, &set);
   err=sched_setaffinity(0, sizeof(set), &set);
   /* Report on rank 0 only without clogging the output */
   if(me==0)
   fprintf(stderr," Core affinity is set to %d on rank %d (not reported on other ranks).\n",cpuid,me);
   /* fprintf(stderr," Core affinity return %d mask %d\n",err,); */
 }

 }

 # endif
