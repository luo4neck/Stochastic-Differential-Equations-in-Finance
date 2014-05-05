#include<iostream>
#include<fstream>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<cmath>
#define N 512 

using namespace std;

const double A[2][2] = {{-1, 1}, {1, -1}};
const double B[2][2] = {{1, 0}, {0, 1}};

void generate_W(gsl_rng *r, double w1[N], double dt)
{
    for(int i=0; i<N; ++i)
    {
       w1[i] = sqrt(dt) * gsl_ran_gaussian(r, 1.0); 
    }
}

int main()
{
gsl_rng *r = gsl_rng_alloc( gsl_rng_mt19937);
gsl_rng_env_setup();
gsl_rng_set(r, time(NULL));

const double T=1.0, dt = T/N;
double dW[N], X[N][2];

ofstream plot("plot.dat");
//for(int j=0; j<MC; ++j) // monte carlo part..
//{
    generate_W(r, dW, dt); // generate two brownian motion paths..
    
    X[0][0] = 1, X[0][1] = 2;
    for(int i=1; i<N; ++i)// approach for euler SDE solution..
    {
        double Xa0 = X[i-1][0]*A[0][0] + X[i-1][1]*A[1][0];
        double Xa1 = X[i-1][0]*A[0][1] + X[i-1][1]*A[1][1];
        
        double Xb0 = X[i-1][0]*B[0][0] + X[i-1][1]*B[1][0];
        double Xb1 = X[i-1][0]*B[0][1] + X[i-1][1]*B[1][1];
        
        X[i][0] = X[i-1][0] + Xa0 * dt + Xb0 * dW[i-1];
        X[i][1] = X[i-1][1] + Xa1 * dt + Xb1 * dW[i-1];
        plot<<i<<" "<<X[i][0]<<" "<<X[i][1]<<endl;
    }
//}
gsl_rng_free(r);
plot.close();


FILE *fp = popen("gnuplot -persist", "w");
if( fp == NULL)
{
    cout<<"Cannot plot the data!"<<endl;
    exit(0);
}
//fprintf(fp, "plot 'plot.dat' u 1:2 w l, 'plot.dat' u 1:3 w l, 'plot.dat' u 1:4 w l, 'plot.dat' u 1:5 w l, 'plot.dat' u 1:6 w l, 'plot.dat' u 1:7 w l, 'plot.dat' u 1:8 w l, 'plot.dat' u 1:9 w l\n");
fprintf(fp, "plot 'plot.dat' using 1:2 w l title 'Euler 0', 'plot.dat' using 1:3 w l title 'Euler 1'\n"); 

return 0;
}
