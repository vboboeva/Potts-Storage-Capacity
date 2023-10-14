#include <stdio.h>
#include <stdlib.h>
#include "random.h"

#define N 100000

int main (){

    int i;
    double x[5], p[5]; /* vectors with values (x) and probability (p) */
    double r[N];

    for (i=0; i<5; i++){
        x[i] = 1.*i;
        p[i] = 1.;
    }

    int seed = rlxd_seed();
    rlxd_init(1,seed);

    discr_dble(r,x,p,5,N);

    for(i=0; i<N; i++) printf("%lf\n",r[i]);

    exit(EXIT_SUCCESS);
}
