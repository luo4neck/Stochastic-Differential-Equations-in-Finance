#include<iostream>
#include<fstream>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<cmath>
#define N 512 

using namespace std;

const double A[2][2] = {{-1, 1}, {1, -1}};
const double B[2][2] = {{1, 0}, {0, 1}};
const double AB[2][2] = {{-1.5, 1}, {1, -1.5}};// AB is the A-0.5*B*B part in the analytical solution.. 

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
    double dW[N], Xeu[N][2], Xan[N][2], Wt=0;
    double ABW[2][2];// ABW will be used in analytical solution..

    ofstream plot("plot.dat");
    generate_W(r, dW, dt); // generate two brownian motion paths..
    
    Xeu[0][0] = 1, Xeu[0][1] = 2, Xan[0][0] = 1, Xan[0][1] = 2;
    Wt = dW[0];
    for(int i=1; i<N; ++i)// approach for euler SDE solution..
    {
        // euler solution..
        double Xa0 = Xeu[i-1][0]*A[0][0] + Xeu[i-1][1]*A[1][0];
        double Xa1 = Xeu[i-1][0]*A[0][1] + Xeu[i-1][1]*A[1][1];
        
        double Xb0 = Xeu[i-1][0]*B[0][0] + Xeu[i-1][1]*B[1][0];
        double Xb1 = Xeu[i-1][0]*B[0][1] + Xeu[i-1][1]*B[1][1];
        
        Xeu[i][0] = Xeu[i-1][0] + Xa0 * dt + Xb0 * dW[i-1];
        Xeu[i][1] = Xeu[i-1][1] + Xa1 * dt + Xb1 * dW[i-1];
    
        //analytic solution..
        Wt = Wt + dW[i];
        
        ABW[0][0] = AB[0][0]*dt*(double)i + B[0][0] * Wt; 
        ABW[0][1] = AB[0][1]*dt*(double)i + B[0][1] * Wt; 
        ABW[1][0] = AB[1][0]*dt*(double)i + B[1][0] * Wt; 
        ABW[1][1] = AB[1][1]*dt*(double)i + B[1][1] * Wt; 
        
        ABW[1][1] = ABW[1][1] - ABW[1][0] * ABW[0][1] / ABW[0][0]; // doing gaussian elimination..
        ABW[1][0] = 0;
        ABW[0][1] = 0;
        //cout<<ABW[0][0]<<" "<<ABW[0][1]<<endl<<ABW[1][0]<<" "<<ABW[1][1]<<endl;

        Xan[i][0] = exp(ABW[0][0]) * Xan[0][0]; 
        Xan[i][1] = exp(ABW[1][1]) * Xan[0][1]; 
        
        plot<<i<<" "<<Xeu[i][0]<<" "<<Xeu[i][1]<<" "<<Xan[i][0]<<" "<<Xan[i][1]<<endl;
    }
    gsl_rng_free(r);
    plot.close();

/* ================= plotting part ================ */

    FILE *fp = popen("gnuplot -persist", "w");
    if( fp == NULL)
    {
        cout<<"Cannot plot the data!"<<endl;
        exit(0);
    }
    fprintf(fp, "set title '2-Dimension SDE solution from Analycital Method and Euler-Maruyama Method'\n");
    fprintf(fp, "set logscale y\n"); 
    fprintf(fp, "plot 'plot.dat' using 1:2 w l title 'Euler 0', 'plot.dat' using 1:3 w l title 'Euler 1', 'plot.dat' using 1:4 w l title 'Actual 0', 'plot.dat' using 1:5 w l title 'Actual 1'\n"); 

    return 0;
}
