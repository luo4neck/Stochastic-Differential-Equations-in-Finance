#include<iostream>
#include<fstream>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<cmath>
#define N 1000

using namespace std;

int main()
{
gsl_rng *R = gsl_rng_alloc( gsl_rng_mt19937);
gsl_rng_env_setup();
gsl_rng_set(R, time(NULL));

int T=1;
double dt = (double)T / (double)N, dW[5][N], W[5][N];
ofstream file("plot.dat");
    for(int i=0; i<N; ++i) {
        for(int j=0; j<5; ++j) {
        if( i == 0 ) 
            { W[j][0] = dW[j][0] = sqrt(dt) * gsl_ran_gaussian(R, 0.5); }
        else{
            dW[j][i] = sqrt(dt) *  gsl_ran_gaussian(R, 0.5);
            W[j][i] = W[j][i-1] + dW[j][i];
            }
        }
    file<<i<<" "<<W[0][i]<<" "<<W[1][i]<<" "<<W[2][i]<<" "<<W[3][i]<<" "<<W[4][i]<<endl;
    }
file.close();
gsl_rng_free(R);

/*=============== plotting part ===============*/
FILE *gp = popen("gnuplot -persist", "w");;

if(gp == NULL) {
printf("Cannot plot the data!\n");
exit(0);
}

fprintf(gp, "set title '5 paths of brownian motion'\n");
fprintf(gp, "plot 'plot.dat' u 1:2 w l title 'path 1', 'plot.dat' u 1:3 w l title 'path 2', 'plot.dat' u 1:4 w l title 'path 3', 'plot.dat' u 1:5 w l title 'path 4', 'plot.dat' u 1:6 w l title 'path 5'\n");
fprintf(gp, "pause -1\n");
fclose(gp);
//system("make clean");
return 0;
}
