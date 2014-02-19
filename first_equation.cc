#include<iostream>
#include<fstream>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<cmath>

using namespace std;

int main()
{
/*gsl_rng *r = gsl_rng_alloc( gsl_rng_mt19937);
gsl_rng_env_setup();
gsl_rng_set(r, time(NULL));

double lambda = 2, mu=1, X[256], dt = 1.0 / (double)N, W[N], Xem[64];
*/
double X[10];
X[0] = 1;

ofstream file("plot.dat");
file<<0<<" "<<1<<endl;
for(int i=1; i<100; ++i)
    {
    X[i] = X[i-1] + (-5.0 * 0.1) * X[i-1]; 
        file<<i<<" "<<X[i]<<endl;
    }
file.close();

/*=============== plotting part ===============*/
FILE *gp = popen("gnuplot -persist", "w");;

if(gp == NULL) {
printf("Cannot plot the data!\n");
exit(0);
}

fprintf(gp, "plot 'plot.dat' u 1:2 w l\n");
fprintf(gp, "pause -1\n");
fclose(gp);
system("cat plot.dat");

//system("make clean");
return 0;
}
