// this file solves dX = aXdt + bXdW with both euler and milstein method.
// a and b are constants, I set a=b=2 here.

#include<iostream>
#include<fstream>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_cdf.h>
#include<cmath>

//#include "function0.h" 
//#include "function1.h" 
#include "function2.h" 
/*
Change header file number to simulate different question in 3.b part.. 
Header file function?.h contains euler function, milstein function and analytic solution fucntion of different SDEsin 3.b part.. 
function0.h is same as 3.a part..
function1.h is 3.b.i part.. 
function2.h is 3.b.ii part.. 
*/

#include "class.h" 
/* class.h
this header file includes two classes and two functions, 
class Line_strong and Line_weak are used to simulate strong simulation and weak simulation..
function double sign() is used to get the negative or positive sign of a double value..
function void Prepare() is used to prepare the random number and true solution of a function..
SIZE 512 and M 50000 are defined in class.h and function header files.. 
*/

using namespace std;

int main()
{
gsl_rng *r = gsl_rng_alloc( gsl_rng_mt19937 );
gsl_rng_env_setup();
gsl_rng_set(r, time(NULL));
ofstream error("error.dat"); 
Line_strong mil3(SIZE), elr3(SIZE), mil2(SIZE/2), elr2(SIZE/2), mil1(SIZE/4), elr1(SIZE/4), mil0(SIZE/8), elr0(SIZE/8);
Line_weak weak3(SIZE), weak2(SIZE/2), weak1(SIZE/4), weak0(SIZE/8);

double sum=0, dt = 1.0/(double) SIZE, Xtrue[SIZE], W[SIZE], dW[SIZE];
int I;

for(int j=0; j<M; ++j) // pathwise loop, M is the time of monte carlo tests..
{

    Prepare(r, dt, dW, W, Xtrue);    
    for(int i=1; i<SIZE; ++i) // stepwise approach
    { // step size 1/512..
        mil3.X[i] = milstein(mil3.X[i-1], dW[i-1], dt); 
        mil3.err[j] = mil3.err[j] + abs(Xtrue[i] - mil3.X[i]);
        
        elr3.X[i] = euler(elr3.X[i-1], dW[i-1], dt);
        elr3.err[j] = elr3.err[j] + abs(Xtrue[i] - elr3.X[i]); 
        
        weak3.X[i]= weak3.X[i-1] + a*weak3.X[i-1]*dt + b*weak3.X[i-1] * dW[i-1];
        //weak3.X[i]= weak3.X[i-1] + a*weak3.X[i-1]*dt + b*weak3.X[i-1]*sign(dW[i-1]) * sqrt(dt); 
        
        if(i%2==0)
        {// step size 1/256..
            I = i/2;
            sum = dW[i-1] + dW[i-2];
            mil2.X[I] = milstein(mil2.X[I-1], sum, dt*2); 
            mil2.err[j] = mil2.err[j] + abs(Xtrue[I*2] - mil2.X[I]);
            
            elr2.X[I] = euler(elr2.X[I-1], sum, 2*dt);
            elr2.err[j] = elr2.err[j] + abs(Xtrue[I*2] - elr2.X[I]);
            
            weak2.X[I]= weak2.X[I-1] + a*weak2.X[I-1]*2*dt + b*weak2.X[I-1] * sum; 
            //weak2.X[I]= weak2.X[I-1] + a*weak2.X[I-1]*2*dt + b*weak2.X[I-1] * sign(sum) * sqrt(2*dt); 
        
            if(i%4==0)
            {// step size 1/128..
                I = i/4;
                sum = dW[i-1] + dW[i-2] + dW[i-3] + dW[i-4];
                mil1.X[I] = milstein(mil1.X[I-1], sum, 4*dt); 
                mil1.err[j] = mil1.err[j] + abs(Xtrue[I*4] - mil1.X[I]);
                
                elr1.X[I] = euler(elr1.X[I-1], sum, 4*dt);
                elr1.err[j] = elr1.err[j] + abs(Xtrue[I*4] - elr1.X[I]);
                
                weak1.X[I]= weak1.X[I-1] + a*weak1.X[I-1]*4*dt + b*weak1.X[I-1] * sum; 
                //weak1.X[I]= weak1.X[I-1] + a*weak1.X[I-1]*4*dt + b*weak1.X[I-1]*sign(sum) * sqrt(4*dt); 
            
                if(i%8==0)
                {// step size 1/64..
                    I = i/8;
                    sum = dW[i-1] + dW[i-2] + dW[i-3] + dW[i-4] + dW[i-5] + dW[i-6] + dW[i-7] + dW[i-8];
                    mil0.X[I] = milstein(mil0.X[I-1], sum, 8*dt); 
                    mil0.err[j] = mil0.err[j] + abs(Xtrue[I*8] - mil0.X[I]);
                    
                    elr0.X[I] = euler(elr0.X[I-1], sum, 8*dt);
                    elr0.err[j] = elr0.err[j] + abs(Xtrue[I*8] - elr0.X[I]);
                    
                    weak0.X[I]= weak0.X[I-1] + a*weak0.X[I-1]*8*dt + b*weak0.X[I-1] * sum; 
                    //weak0.X[I]= weak0.X[I-1] + a*weak0.X[I-1]*8*dt + b*weak0.X[I-1]*sign(sum) * sqrt(8*dt); 
                }
            }
        }
    }
mil3.err[j] = mil3.err[j]/(double)mil3.size;
mil2.err[j] = mil2.err[j]/(double)mil2.size;
mil1.err[j] = mil1.err[j]/(double)mil1.size;
mil0.err[j] = mil0.err[j]/(double)mil0.size;
mil3.mean = mil3.mean + mil3.err[j]; //accumulate the mean..
mil2.mean = mil2.mean + mil2.err[j]; 
mil1.mean = mil1.mean + mil1.err[j];
mil0.mean = mil0.mean + mil0.err[j];

elr3.err[j] = elr3.err[j]/(double)elr3.size;
elr2.err[j] = elr2.err[j]/(double)elr2.size;
elr1.err[j] = elr1.err[j]/(double)elr1.size;
elr0.err[j] = elr0.err[j]/(double)elr0.size;
elr3.mean = elr3.mean + elr3.err[j];
elr2.mean = elr2.mean + elr2.err[j];
elr1.mean = elr1.mean + elr1.err[j];
elr0.mean = elr0.mean + elr0.err[j];

weak3.mean = weak3.mean + weak3.X[weak3.size-1];
weak2.mean = weak2.mean + weak2.X[weak2.size-1];
weak1.mean = weak1.mean + weak1.X[weak1.size-1];
weak0.mean = weak0.mean + weak0.X[weak0.size-1];

error<<j<<" "<<mil3.err[j]<<" "<<elr3.err[j]<<" "<<mil2.err[j]<<" "<<elr2.err[j]<<" "<<mil1.err[j]<<" "<<elr1.err[j]<<" "<<mil0.err[j]<<" "<<elr0.err[j]<<endl;// print the error of two methods via gnuplot..

}// loop over..

gsl_rng_free(r);
error.close();

mil3.mean = mil3.mean/(double)M;//get the meam..
mil2.mean = mil2.mean/(double)M;
mil1.mean = mil1.mean/(double)M;
mil0.mean = mil0.mean/(double)M;
elr3.mean = elr3.mean/(double)M;
elr2.mean = elr2.mean/(double)M;
elr1.mean = elr1.mean/(double)M;
elr0.mean = elr0.mean/(double)M;

/* ===========  weak convergence computing ==================== */

weak3.mean = abs(weak3.mean/(double)M - exp(a));// weak convergence..
weak2.mean = abs(weak2.mean/(double)M - exp(a));
weak1.mean = abs(weak1.mean/(double)M - exp(a));
weak0.mean = abs(weak0.mean/(double)M - exp(a));

/* ===================== weak part ====================== */

ofstream data("data.dat");// send all data into 'data.dat' for error convergence..
data<<dt  <<" "<<dt  <<" "<<mil3.mean<<" "<<mil3.up(0.975)<<" "<<mil3.down(0.975)<<" "<<mil3.up(0.95)<<" "<<mil3.down(0.95)<<" "<<elr3.mean<<" "<<elr3.up(0.975)<<" "<<elr3.down(0.975)<<" "<<elr3.up(0.95)<<" "<<elr3.down(0.95)<<" "<<weak3.mean<<endl;

data<<dt*2<<" "<<dt*2<<" "<<mil2.mean<<" "<<mil2.up(0.975)<<" "<<mil2.down(0.975)<<" "<<mil2.up(0.95)<<" "<<mil2.down(0.95)<<" "<<elr2.mean<<" "<<elr2.up(0.975)<<" "<<elr2.down(0.975)<<" "<<elr2.up(0.95)<<" "<<elr2.down(0.95)<<" "<<weak2.mean<<endl;

data<<dt*4<<" "<<dt*4<<" "<<mil1.mean<<" "<<mil1.up(0.975)<<" "<<mil1.down(0.975)<<" "<<mil1.up(0.95)<<" "<<mil1.down(0.95)<<" "<<elr1.mean<<" "<<elr1.up(0.975)<<" "<<elr1.down(0.975)<<" "<<elr1.up(0.95)<<" "<<elr1.down(0.95)<<" "<<weak1.mean<<endl;

data<<dt*8<<" "<<dt*8<<" "<<mil0.mean<<" "<<mil0.up(0.975)<<" "<<mil0.down(0.975)<<" "<<mil0.up(0.95)<<" "<<mil0.down(0.95)<<" "<<elr0.mean<<" "<<elr0.up(0.975)<<" "<<elr0.down(0.975)<<" "<<elr0.up(0.95)<<" "<<elr0.down(0.95)<<" "<<weak0.mean<<endl;
data.close();

/* ========== data store part ended, error plot part start ========== */

FILE *fp = popen("gnuplot -persist", "w");
if(fp == NULL)
{
    cout<<"Cannot plot the data!"<<endl;
    exit(0);
}
fprintf(fp, "set title '8 Line show the error of 4 step sizes and 2 method'\n");
fprintf(fp, "unset key\n");
fprintf(fp, "set xlabel 'Number of batches'\n");
fprintf(fp, "set ylabel 'Errors at different batches'\n");
fprintf(fp, "plot 'error.dat' u 1:2 w l, 'error.dat' u 1:3 w l, 'error.dat' u 1:4 w l, 'error.dat' u 1:5 w l, 'error.dat' u 1:6 w l, 'error.dat' u 1:7 w l, 'error.dat' u 1:8 w l, 'error.dat' u 1:9 w l\n");
fprintf(fp, "pause -1\n");
fclose(fp);

/* ========== error plot ended, error bar plot start ========== */

fp = popen("gnuplot -persist", "w");
if(fp == NULL)
{
    cout<<"Cannot plot the data!"<<endl;
    exit(0);
}
fprintf(fp, "set title 'Errorbars at different step sizes for Milstein method and Euler-Maruyama method'\n");
fprintf(fp, "set logscale xy\n");
fprintf(fp, "set xlabel 'Step sizes: of 1/512 1/256 1/128 1/64'\n");
fprintf(fp, "set ylabel 'Errorbars at different step sizes'\n");
fprintf(fp, "set xrange [0.001:0.05]\n");
fprintf(fp, "plot 'data.dat' u 1:2 w l title 'standard line', 'data.dat' u 1:3:4:5 w yerrorlines title 'milstein 95 confidence interval', 'data.dat' u 1:3:6:7 w yerrorlines title 'milstein 90 confidence interval', 'data.dat' u 1:8:9:10 w yerrorlines title 'euler 95 confidence interval', 'data.dat' u 1:8:11:12 w yerrorlines title 'euler 90 confidence interval', 'data.dat' u 1:13 w l title 'euler weak convergence'\n");
fprintf(fp, "pause -1\n");
fclose(fp); 

/* ========== error plot ended, fitting line plot start ========== */

fp = popen("gnuplot -persist", "w");
if(fp == NULL)
{
    cout<<"Cannot plot the data!"<<endl;
    exit(0);
}
fprintf(fp, "set title 'Fitting Line and Error at different step sizes for Milstein method and Euler-Maruyama method'\n");
fprintf(fp, "set logscale xy\n");
fprintf(fp, "set xlabel 'Step sizes: of 1/512 1/256 1/128 1/64'\n");
fprintf(fp, "set ylabel 'Errors at different step sizes'\n");
fprintf(fp, "set xrange [0.001:0.05]\n");
fprintf(fp, "f1(x)=a1*x**b1\n");//milstein line..
fprintf(fp, "a1=3;b1=1\n");
fprintf(fp, "f2(x)=a2*x**b2\n");//euler line..
fprintf(fp, "a2=5;b2=1\n");
fprintf(fp, "f3(x)=a3*x**b3\n");//euler weak line..
fprintf(fp, "a3=5;b3=1\n");
fprintf(fp, "fit [0.001:0.05] f1(x) 'data.dat' u 1:3 via a1,b1\n");
fprintf(fp, "fit [0.001:0.05] f2(x) 'data.dat' u 1:8 via a2,b2\n");
fprintf(fp, "fit [0.001:0.05] f3(x) 'data.dat' u 1:13 via a3,b3\n");
fprintf(fp, "plot 'data.dat' u 1:2 w l title 'standard line', 'data.dat' u 1:3 w l title 'error of milstein', 'data.dat' u 1:8 w l title 'error of euler strong', 'data.dat' u 1:13 w l title 'error of euler weak', f1(x) lw 1 lc rgb 'orange' title 'fitting of milstein strong', f2(x) lw 1 lc rgb 'yellow' title 'fitting of euler strong', f3(x) lw 1 lc rgb 'grey' title 'fitting of euler weak'\n");

fprintf(fp, "pause -1\n");
fclose(fp); 
return 0;
}
