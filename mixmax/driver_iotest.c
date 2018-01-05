#include <stdio.h>
#include <stdlib.h>
#include "mixmax.h"

int main( int argc,  char *argv[] ){
    rng_state_t* X=rng_alloc();
    rng_state_t* Y=rng_alloc();
    
    seed_uniquestream(X, 0,  0,  0,  1);

    FILE* fout;
    if (! (fout = fopen("seed1.conf", "w"))){printf("trouble opening file for writing the state\n");};
    X->fh = fout;
    print_state(X);
    fclose(X->fh );
    read_state(Y, "seed1.conf");
    printf("Successfully read the state from seed1.conf\n");

    int i;
    for(i=0;i<N;i++){
        if (X->V[i] != Y->V[i]) {
            printf("Reading Error: component %d does not match %llu !=%llu\n", i, X->V[i], Y->V[i]);
            exit(-1);
        }
    }
    printf("State Reading Test ok\n");

    rng_free(X);
    rng_free(Y);
	return 0;
}
