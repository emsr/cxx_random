/*
 
  MIXMAX generator - Example usage
 
 */

#include <stdio.h>
#include <stdlib.h>
#include "mixmax.h"

int main( int argc,  char *argv[] ){
	myuint j,p;
	//double array[ARRAY_SIZE];
	
    fprintf(stderr,"Welcome to the MIXMAX random number generator!\nThe curent matrix size is N=%u\n"
            "(the actual matrix is not kept in memory in this new efficient implementation)\n"
            "special entry in the matrix is %lld\n"
            "special multiplier m=2^%u+1\n"
            "Working in the Galois field with modulus 2^%u-1\n", rng_get_N(), (int64_t)SPECIAL, SPECIALMUL, BITS);
	
	fprintf(stderr,"\nHow many numbers?\n"
                   "(2^10=1024, 2^20=1048576, 2^30=1073741824)\nEnter m: ");
	if (!scanf("%llu", &p)) {printf("error reading number"); exit(-1);}
	
	rng_state_t S;
	rng_state_t* X = &S;
	// or allocate dynamically:
	//rng_state_t* X=rng_alloc();
	
	//seed_spbox(X, 123);      // seed with nonlinear SP-box, guarantees a unique seed and independence, but not non-collision of different streams, seed can be 1 ... 2^64-1
	

	/* Best seeding method here - guarantees complete independence,
       uniqueness and non-colision of different streams!
	   just do:                                                          */
    // seed_uniquestream(X, clusterID,  machineID,  runID,   streamID);
    seed_uniquestream(X, 0,  0,  0,  1);
	
	for (j=0; j<p ; j++) {
        printf("%1.18F\n", get_next_float(X) );
//               printf("%1.18F %1.18F %1.18F %1.18F %1.18F\n", get_next_float(X), get_next_float(X),
//                      get_next_float(X), get_next_float(X), get_next_float(X) );
        // for floating point number on (0,1]
	}
    fprintf(stdout, "ok\n");
    
	// if state was allocated dynamically, free it when its no longer needed:
	//rng_free(S);
	return 0;
}
