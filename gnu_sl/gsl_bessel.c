#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_math.h>


int main(){
    double x = 5.0;
    double y = gds_sf_bessel_J0(x);
    printf ("J0(%g) = %0.18e\n", x, y);

    return 0;
}
