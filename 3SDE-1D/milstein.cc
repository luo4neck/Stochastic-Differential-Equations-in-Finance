// this file solves dX = aXdt + bXdW with milstein method.
// a and b are constants, I set a=b=2 here.

// X0 is analytic solution and X1, X2, X3 are numercial solution..
// X0 is red line with dt = 0.001
// X1 is green line with dt = 0.001
// X2 is blue line with dt = 0.01
// X3 is purple line with dt = 0.1

#include<iostream>
#include<fstream>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<cmath>
#define MAX 1000

using namespace std;

int main()
{
gsl_rng *r = gsl_rng_alloc( gsl_rng_mt19937 );
gsl_rng_env_setup();
gsl_rng_set(r, time(NULL));

double X1[MAX], a=2, b=2, dt = 1.0/(double) MAX, W[MAX], dW[MAX];
//double X0[MAX], X2[MAX/10], X3[MAX/100]

X1[0] = 1;
W[0] = dW[0] = sqrt(dt) * gsl_ran_gaussian(r, 1.0);

ofstream file("plot.dat");
    for(int i=1; i<MAX; ++i)
    {
        dW[i] = sqrt(dt) * gsl_ran_gaussian(r, 1.0);
        W[i] = W[i-1] + dW[i-1];
        X1[i]= X1[i-1] + a*X1[i-1]*dt + b*X1[i-1]*dW[i-1] + 0.5*b*X1[i-1]*b*( dW[i-1]*dW[i-1] - dt);
        if((i+1)%10 == 0)
        {
            if((i+1)%100 == 0)
            {
            file<<(double)i/MAX<<" "<<X1[i]<<" "<<X1[i]<<" "<<X1[i]<<" "<<0<<endl;
            }
            else 
            file<<(double)i/MAX<<" "<<X1[i]<<" "<<X1[i]<<endl;

        }
        else 
        file<<(double)i/MAX<<" "<<X1[i]<<endl;
    }
gsl_rng_free(r);
file.close();

/* ========== data store part ended, plotting part start ========== */

FILE *gp = popen("gnuplot -persist", "w");

if(gp == NULL)
{
cout<<"Cannot plot the data!"<<endl;
exit(0);
}

//fprintf(gp, "set logscale xy\n");
fprintf(gp, "plot 'plot.dat' u 1:2 w l, 'plot.dat' u 1:3 w l, 'plot.dat' u 1:4 w l, 'plot.dat' u 1:5 w l\n");
fprintf(gp, "pause -1\n");
fclose(gp);

return 0;
}
