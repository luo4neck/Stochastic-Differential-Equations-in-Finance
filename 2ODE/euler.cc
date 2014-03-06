#include<iostream>
#include<fstream>
#include<gsl/gsl_rng.h>
#include<cmath>

using namespace std;

//X0 is analyst solution and X1, X2, X3 are euler solution

//X0 is red line is count by i, dt=0.001
//X1 is green line is count by i, dt=0.001
//X2 is blue line is count by j, dt=0.01
//X3 is purple line is count by k, dt=0.1
int main()
{
int j=1, k=1;
double X0[1000],X1[1000], X2[100], X3[10];
X1[0] = X0[0] = X2[0] =  X3[0] = 1;

ofstream file("plot.dat");
file<<0<<" "<<1<<" "<<1<<" "<<1<<" "<<1<<endl;
for(int i=1; i<1000; ++i)
    {
    X0[i] = exp( -5 * 0.001 * i);
    X1[i] = X1[i-1] - 5.0 * 0.001 * X1[i-1]; 
    if( (i)%10 == 0)
        {
        X2[j] = X2[j-1] - 5.0 * 0.01  * X2[j-1];
        if ( (i)%100 == 0)
            {
            X3[k] = X3[k-1] - 5.0 * 0.1  * X3[k-1];
            file<<i<<" "<<X0[i]<<" "<<X1[i]<<" "<<X2[j]<<" "<<X3[k]<<endl;
            k++;
            }
        else{
            file<<i<<" "<<X0[i]<<" "<<X1[i]<<" "<<X2[j]<<endl;
            }
        j++;
        }
    else{
        file<<i<<" "<<X0[i]<<" "<<X1[i]<<endl;
        }
    }
file.close();

/*=============== plotting part ===============*/
FILE *gp = popen("gnuplot -persist", "w");;

if(gp == NULL) {
printf("Cannot plot the data!\n");
exit(0);
}

fprintf(gp, "plot 'plot.dat' u 1:2 w l, 'plot.dat' u 1:3 w l, 'plot.dat' u 1:4 w l, 'plot.dat' u 1:5 w l\n");
fprintf(gp, "pause -1\n");
fclose(gp);
system("vim plot.dat");
return 0;
}
