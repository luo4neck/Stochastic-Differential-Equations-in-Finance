#include<iostream>
#include<fstream>
#include<gsl/gsl_rng.h>
#include<cmath>

using namespace std;

int main()
{
int j=1, k=1;//j is for dt=0.01 and X2, k is for dt=0.1 and X3
double X0[1000],X1[1000], X2[100], X3[10];
X1[0] = X0[0] = X2[0] =  X3[0] = 1;

ofstream file("plot.dat");
file<<0<<" "<<1<<" "<<1<<" "<<1<<" "<<1<<endl;
for(int i=1; i<1000; ++i)
    {
    if( (i)%10 == 0)
        {
        if ( (i)%100 == 0)
            {
            X3[k] = X3[k-1] - 5.0 * 0.1  * X3[k-1];
            X2[j] = X2[j-1] - 5.0 * 0.01  * X2[j-1];
            X1[i] = X1[i-1] - 5.0 * 0.001 * X1[i-1]; 
            X0[i] = exp( -5 * 0.001 * i);
            file<<i<<" "<<X0[i]<<" "<<X1[i]<<" "<<X2[j]<<" "<<X3[k]<<endl;
            k++, j++;
            }
        else{
            X2[j] = X2[j-1] - 5.0 * 0.01  * X2[j-1];
            X1[i] = X1[i-1] - 5.0 * 0.001 * X1[i-1]; 
            X0[i] = exp( -5 * 0.001 * i);
            file<<i<<" "<<X0[i]<<" "<<X1[i]<<" "<<X2[j]<<endl;
            j++;
            }
        }
    else{
        X1[i] = X1[i-1] - 5.0 * 0.001 * X1[i-1]; 
        X0[i] = exp( -5 * 0.001 * i);
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
return 0;
}
