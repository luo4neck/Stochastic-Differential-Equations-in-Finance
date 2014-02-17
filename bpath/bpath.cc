#include<iostream>
#include<fstream>
#include<gsl/gsl_rng.h>

using namespace std;

int main()
{
const gsl_rng* R;
gsl_rng_env_setup();
R = gsl_rng_alloc( gsl_rng_mt19937);
gsl_rng_set(R, time(NULL));

int T=1, N = 500;
double dt = (double)T / (double)N, dW[500], W[500];

for(int i=0; i<500; ++i) {
dW[i] = W[i] = 0;
}

ofstream file("plot.dat");
    for(int i=0; i<100; ++i) {
    double a = gsl_rng_uniform(R);
    cout<<a<<", ";
    if ( (i+1)%10 == 0) cout<<endl;
    file<<i<<" "<<a<<endl;
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
system("make clean");
return 0;
}
