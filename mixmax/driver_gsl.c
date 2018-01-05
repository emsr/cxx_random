/*
 Example usage of MIXMAX with GSL
 
 */

// cc -std=c99 -O3 -funroll-loops -Wall -I /opt/local/include/ driver_gsl.c mixmax.o -L/opt/local/lib -DHOOKUP_GSL=1  -o gsl -lgsl
// cc -std=c99 -O3 -funroll-loops -Wall -I /usr/local/include/ driver_gsl.c mixmax.o -L/usr/local/lib -DHOOKUP_GSL=1  -o gsl -lgsl
// make gsl; GSL_RNG_SEED=123 GSL_RNG_TYPE=MIXMAX ./gsl

#include <stdio.h>
#include "mixmax.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_version.h>

int
main (void)
{
  const gsl_rng_type * T;
  gsl_rng * r;
  unsigned long int seed = 0;

  int i, n = 10;

 // gsl_rng_env_setup();
 // gsl_rng_types_setup ();

  T = gsl_rng_mixmax;
    
    const char *p = getenv ("GSL_RNG_SEED");
    if (p)
    {
        seed = strtoul (p, 0, 0);
        fprintf (stderr, "GSL_RNG_SEED=%lu\n", seed);
    };
    
    gsl_rng_default_seed = seed;

    r = gsl_rng_alloc (T);

    printf("generator type: %s\n", gsl_rng_name (r));
    printf("seed = %lu\n", gsl_rng_default_seed);
    printf("gsl_version=%s\n", gsl_version);
    
  for (i = 0; i < n; i++) 
    {
      double u = gsl_rng_uniform (r);
      printf ("%.18F\n", u);
    }

  gsl_rng_free (r);

  return 0;
}
