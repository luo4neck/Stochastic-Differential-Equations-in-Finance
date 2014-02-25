#include<iostream>
#include<fstream>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<cmath>
#define N 512 
#define M 1000
using namespace std;

int main()
{
gsl_rng *r = gsl_rng_alloc( gsl_rng_mt19937);
gsl_rng_env_setup();
gsl_rng_set(r, time(NULL));

double lambda = 2, mu=1, Xerr[5][M], T=1, dt = 1.0 / (double)N;
double W[M], dW[M], Xtrue[M];

for(int i=0; i<M; ++i)
    {   
    for(int j=0; j<5; ++j) Xerr[j][i] = 0;
    }

ofstream file("plot.dat");
for(int i=0; i<M; ++i)
    {
    if(i==0)
        {
        W[i] = dW[i] = gsl_ran_gaussian(r, 1.0);
        Xtrue[i] = 1;
        }
    else{
        dW[i] = sqrt(dt) * gsl_ran_gaussian(r, 1.0);
        W[i] = W[i-1] + dW[i-1];
        Xtrue[i] = Xtrue[0] * exp( (lambda-0.5*mu*mu)*dt*i + mu*W[i] );
        }

    for(int p=1; p<=5; ++p)
        {
        double R = pow(2, p-1), Dt = dt * R, L = N/R, Xtmp = Xtrue[0];
        for(int j=1; j<=L; ++j)
            {
                double Winc = 0;
   //             cout<<i<<","<<p<<","<<j<<" "<<"well here!"<<endl;
                for(int k=R*(j-1)+1; k<=R*(j-1)+R*j; ++k)
                Winc = Winc + dW[k];
                Xtmp = Xtmp + Dt * lambda * Xtmp + mu* Xtmp * Winc;
            }
     //           cout<<"oh no!"<<endl;
        Xerr[p-1][i] = abs(Xtmp - Xtrue[i]) ;
        }
    file<<i<<" "<<Xerr[0][i]<<" "<<Xerr[1][i]<<" "<<Xerr[2][i]<<" "<<Xerr[3][i]<<" "<<Xerr[4][i]<<endl;
    }
gsl_rng_free(r);
file.close();

/*
for(int i=0; i<M; ++i)
    {   
    cout<<Xerr[0][i]<<" "<<Xerr[1][i]<<" "<<Xerr[2][i]<<" "<<Xerr[3][i]<<" "<<Xerr[4][i]<<endl;
    }

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
*/
/*=============== plotting part ===============*/


FILE *gp = popen("gnuplot -persist", "w");;

if(gp == NULL) {
printf("Cannot plot the data!\n");
exit(0);
}

fprintf(gp, "plot 'plot.dat' u 1:2 w l, 'plot.dat' u 1:3 w l, 'plot.dat' u 1:4 w l, 'plot.dat' u 1:5 w l, 'plot.dat' u 1:6 w l\n");
fprintf(gp, "pause -1\n");
fclose(gp);
//system("cat plot.dat");
//system("make clean");

return 0;
}
