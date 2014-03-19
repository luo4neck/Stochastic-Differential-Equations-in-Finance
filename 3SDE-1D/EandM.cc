// this file solves dX = aXdt + bXdW with both euler and milstein method.
// a and b are constants, I set a=b=2 here.

// X0 is red line with dt = 0.001 | milstein method
// X1 is green line with dt = 0.001 | euler method

// X2 is blue line with dt = 0.01 | milstein method
// X3 is purple line with dt = 0.01 | euler method

// X4 is purple line with dt = 0.1 | milstein method
// X5 is purple line with dt = 0.1 | euler method

#include<iostream>
#include<fstream>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_cdf.h>
#include<cmath>
#define SIZE 512
#define M 100

using namespace std;

double Prepare(gsl_rng *r, double dt, double a, double b, double dW[SIZE], double W[SIZE], double Xtrue[SIZE])
{
    double sum=Xtrue[0];
    W[0] = dW[0] = sqrt(dt) * gsl_ran_gaussian(r, 0);
    for(int i=1; i<SIZE; ++i)
    {
        dW[i] = sqrt(dt) * gsl_ran_gaussian(r, 1.0);
        W[i] = W[i-1] + dW[i-1];
        Xtrue[i] = 1.0 * exp( (a-0.5*b*b) + b*W[i]);
        sum = sum + Xtrue[i];
    }
    return sum;
}

int main()
{
gsl_rng *r = gsl_rng_alloc( gsl_rng_mt19937 );
gsl_rng_env_setup();
gsl_rng_set(r, time(NULL));
ofstream file("plot.dat");

double X0[SIZE], X1[SIZE], a=2, b=2, dt = 1.0/(double) SIZE, Xtrue[SIZE], W[SIZE], dW[SIZE];
double sumXtrue, err0[M], err1[M];
double mean0=0, mean1=0;
//double X0[SIZE], X2[SIZE/10], X3[SIZE/100]

for(int j=0; j<M; ++j)
{
Xtrue[0] = X0[0] = X1[0] = 1; 
err0[j] = err1[j] = 0;
sumXtrue = Prepare(r, dt, a, b, dW, W, Xtrue);    
    for(int i=1; i<SIZE; ++i)
    {
        X0[i]= X0[i-1] + a*X0[i-1]*dt + b*X0[i-1]*dW[i-1] + 0.5*b*X0[i-1]*b*( dW[i-1]*dW[i-1] - dt);// milstein method..
        err0[j] = err0[j] + abs(Xtrue[i] - X0[i]);
        
        X1[i]= X1[i-1] + a*X1[i-1]*dt + b*X1[i-1]*dW[i-1]; // euler-maruyama method.. 
        err1[j] = err1[j] + abs(Xtrue[i] - X1[i]);
     /*   if((i+1)%10 == 0)
        {
            if((i+1)%100 == 0)
            {
            file<<(double)i/SIZE<<" "<<X1[i]<<" "<<X1[i]<<" "<<X1[i]<<" "<<0<<endl;
            }
            else 
            file<<(double)i/SIZE<<" "<<X1[i]<<" "<<X1[i]<<endl;
        }
        else*/
//file<<(double)i/SIZE<<" "<<X0[i]<<" "<<X1[i]<<" "<<Xtrue[i]<<endl;
    }
err0[j] = err0[j]/(double)SIZE;
err1[j] = err1[j]/(double)SIZE;
mean0 = mean0 + err0[j];
mean1 = mean1 + err1[j];
file<<j<<" "<<err0[j]<<" "<<err1[j]<<endl;
//cout<<sumXtrue<<endl<<sumX0<<endl<<sumX1<<endl;
}
gsl_rng_free(r);
file.close();
mean0 = mean0/(double)M;
mean1 = mean1/(double)M;

cout<<mean0<<endl<<mean1<<endl;
double sigma0=0, sigma1=0;
for(int i=0; i<M; ++i)
{
    sigma0 = sigma0 + (mean0 - err0[i])*(mean0 - err0[i]);
    sigma1 = sigma1 + (mean1 - err1[i])*(mean1 - err1[i]);
}    
sigma0 = sigma0/(double)(M-1);
sigma1 = sigma1/(double)(M-1);
cout<<sigma0<<endl<<sigma1<<endl;

cout<<"The 95% confidence interval of milstein is:"<<endl;
cout<<mean0 + gsl_cdf_tdist_Pinv(0.95, M-1)*sqrt(sigma0/(double)M)<<" : "<<mean0 + gsl_cdf_tdist_Qinv(0.95, M-1)*sqrt(sigma0/(double)M)<<endl;

cout<<"The 95% confidence interval of euler is:"<<endl;
cout<<mean1 + gsl_cdf_tdist_Pinv(0.95, M-1)*sqrt(sigma1/(double)M)<<" : "<<mean1 + gsl_cdf_tdist_Qinv(0.95, M-1)*sqrt(sigma1/(double)M)<<endl;
/* ========== data store part ended, plotting part start ========== */

FILE *gp = popen("gnuplot -persist", "w");

if(gp == NULL)
{
cout<<"Cannot plot the data!"<<endl;
exit(0);
}

//fprintf(gp, "set logscale xy\n");
//fprintf(gp, "plot 'plot.dat' u 1:2 w l, 'plot.dat' u 1:3 w l, 'plot.dat' u 1:4 w l, 'plot.dat' u 1:5 w l\n");
//fprintf(gp, "plot 'plot.dat' u 1:2 w l, 'plot.dat' u 1:3 w l, 'plot.dat' u 1:4 w l\n");
fprintf(gp, "plot 'plot.dat' u 1:2 w l, 'plot.dat' u 1:3 w l\n");
fprintf(gp, "pause -1\n");
fclose(gp);

//cout<<gsl_cdf_tdist_Pinv(0.95, 19)<<endl; // 0.95 here indicate one side present value and 19 is sample number-1;
//cout<<gsl_cdf_tdist_Qinv(0.95, 19)<<endl; // 0.95 here indicate one side present value and 19 is sample number-1;
//one side 90 95 97.5
//two side 80 90 95
return 0;
}
