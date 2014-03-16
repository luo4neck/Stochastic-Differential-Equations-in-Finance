#include<iostream>
#include<fstream>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<cmath>
//#define N 512 
#define M 50000
using namespace std;

int main()
{
gsl_rng *r = gsl_rng_alloc( gsl_rng_mt19937);
gsl_rng_env_setup();
gsl_rng_set(r, time(NULL));

double lambda = 2, mu=0.1,  T=1, X0=1, Xem[5]={0, 0, 0, 0, 0}, Xtmp[M];
for(int p=1; p<=5; ++p)
    {
    double Dt=pow(2.0, p-10); 
    double L=T/Dt;
    
    for(int i=0; i<M; ++i)
        {Xtmp[i] = X0*1; }
    
    for(int j=1; j<=L; ++j)
        {
        double Winc[M]; 
        for(int i=0; i<M; ++i)
            {
            Winc[i] = sqrt(Dt) * gsl_ran_gaussian(r, 0);
            }
        
        for(int i=0; i<M; ++i)
            {
            Xtmp[i] = Xtmp[i] + Dt * lambda * Xtmp[i] + mu * Xtmp[i] * Winc[i];
            }
        }

     double sum=0;
     for(int i=0; i<M; ++i)
        {
        sum = sum + Xtmp[i];
        }
    Xem[p-1] = sum/(double)M;
    }
gsl_rng_free(r);

double Dtvals[5], Xerr[5];
ofstream file("plot.dat");
for(int i=0; i<5; ++i)
    {
    Xerr[i] = abs(Xem[i] - exp(lambda));
    cout<<Xerr[i]<<endl;
    Dtvals[i] = pow(2, i-9);
    file<<Dtvals[i]<<" "<<Xerr[i]<<" "<<Dtvals[i]<<endl;
    }
file.close();

/*=============== plotting part ===============*/


FILE *gp = popen("gnuplot -persist", "w");;

if(gp == NULL) 
{
    printf("Cannot plot the data!\n");
    exit(0);
}

fprintf(gp, "set logscale xy\n");
fprintf(gp, "plot 'plot.dat' u 1:2 w l, 'plot.dat' u 1:3 w l\n");
fprintf(gp, "pause -1\n");
fclose(gp);
//system("cat plot.dat");
//system("make clean");

return 0;
}
