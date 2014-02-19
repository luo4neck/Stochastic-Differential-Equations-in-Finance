#include<iostream>
#include<fstream>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<cmath>
#define N 256 

using namespace std;

int main()
{
gsl_rng *r = gsl_rng_alloc( gsl_rng_mt19937);
gsl_rng_env_setup();
gsl_rng_set(r, time(NULL));

double lambda = 2, mu=1, X[256], dt = 1.0 / (double)N, W[N], Xem[64];
ofstream file("plot.dat");
W[0] = gsl_ran_gaussian(r, 1.0);
X[0] = 1;
file<<0<<" "<<X[0]<<endl;

int j=0;
for(int i=1; i<N; ++i)
    {
    W[i] = W[i-1] + sqrt(dt) * gsl_ran_gaussian(r, 1.0); 
    X[i] = X[0] * exp( (lambda-0.5*mu*mu)*dt*i + mu*W[i] );
    if( (i+1)%4 == 0)
        {
        Xem[j] = X[i];
        cout<<X[i]<<" "<<Xem[j]<<endl;
        file<<i*dt<<" "<<X[i]<<" "<<Xem[j]<<endl;
        j++;
        }
    else 
        file<<i*dt<<" "<<X[i]<<endl;
    }
file.close();
gsl_rng_free(r);

/*=============== plotting part ===============*/
FILE *gp = popen("gnuplot -persist", "w");;

if(gp == NULL) {
printf("Cannot plot the data!\n");
exit(0);
}

fprintf(gp, "plot 'plot.dat' u 1:2 w l, 'plot.dat' u 1:3 w l\n");
fprintf(gp, "pause -1\n");
fclose(gp);
//system("cat plot.dat");
system("make clean");
return 0;
}
