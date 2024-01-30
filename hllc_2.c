#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

int x;
int N = 1002;
double y = (double) 5/3;
double t, dx;

/* functions to find the minimum and maximum values of an array */

double max(double array[], int n)
{
    int i;
    double max = fabs(array[0]);
    for (i = 1; i < n; i++)
        if (fabs(array[i]) > max)
            max = fabs(array[i]);
    return max;
}

double min(double array[], int n)
{
    int i;
    double min = fabs(array[0]);
    for (i = 1; i < n; i++)
        if (fabs(array[i]) < min)
            min = fabs(array[i]);
    return min;
}

/* functions to find the fluid states and place the correct states into the correct arrays */

void q(double d, double v, double e, double k[])
{
    k[0] = d;
    k[1] = d*v;
    k[2] = e;
}

/* functions to find the flux states and place the correct states into the correct arrays */

void f(double d, double v, double e, double p, double k[])
{
    k[0] = d*v;
    k[1] = d*v*v + p;
    k[2] = v*(e + p);
}

/* function to calcuklate the final fluid states */

double qhfinal(double dt, double fr, double fl, double qm, double dx)
{
    double q = qm - (dt/dx)*(fr - fl);
    return q;
}

/* function to estimate the pressure of the region in between the waves*/

double pressure_estimate(double pl, double pr, double vl, double vr, double dl, double dr, double cl, double cr){
    double pressure_est = fmax(0, 0.5*(pl + pr) - 0.5*(vr - vl)*0.5*(dr + dl)*0.5*(cl + cr));
    return pressure_est;
}

/* function to estimate the sound speed */

double cs(double p, double d){
    double cs = sqrt(y*p*(1/d));
    return cs;
}

/* function to calculate the total energy of the system */

double total_e(double e, double d, double v)
{
    double total_e = e/d - 0.5*v*v;
    return total_e;
}

/* part of the shock wave speed calculation */

double shock_wave(double cs, double pcd, double p){
    double shock = cs*sqrt(1+((1+y)/(2*y))*((pcd/p)-1));
    return shock;
}

/* function to calculate the first intermediate fluid state */

double qmid1(double dl, double vl, double sl, double smid)
{
   double q1 = dl*(sl - vl)/(sl - smid);
   return q1;
}

/* function to calculate the second intermediate fluid state */

double qmid2(double dl, double vl, double sl, double smid)
{
   double q2 = smid*dl*(sl - vl)/(sl - smid);
   return q2;
}

/* function to calculate the third intermediate fluid state */

double qmid3(double dl, double vl, double smid, double sl, double el, double pl)
{
   double q3 = (dl*(sl - vl)/(sl - smid))*(el/dl + (smid - vl)*(smid + (pl/(dl*(sl - vl)))));
   return q3;
}

/* Setting out the variables to be used in the solver */

double fh1, fh2, fh3, denl, denvell, enel, int_el, presl, vell, denr, denvelr, ener, int_er, presr, velr, dt, cl, cr, pres, sr, sl, smid;

int main(){
    t = 0;
    
    /* creating a file and the file pointer */
    FILE *data;
    data = fopen("hllc_2.txt", "w");
    /* allocating the memory for the arrays */
    double *p = (double*) malloc(N * sizeof(double));
    double *d = (double*) malloc(N * sizeof(double));
    double *v = (double*) malloc(N * sizeof(double));
    double *e = (double*) malloc(N * sizeof(double));
    double *energy = (double*) malloc(N * sizeof(double));
    double *q1 = (double*) malloc(N * sizeof(double));
    double *q2 = (double*) malloc(N * sizeof(double));
    double *q3 = (double*) malloc(N * sizeof(double));
    double *q1n = (double*) malloc(N * sizeof(double));
    double *q2n = (double*) malloc(N * sizeof(double));
    double *q3n = (double*) malloc(N * sizeof(double));
    
    double fux[3];
    double flux[3][N];
    double qmid[3][N];
    
    double ql[3];
    double qr[3];
    
    double qlm[3];
    double qrm[3];
    
    double fl[3];
    double fr[3];

    double dx = 0.01;

    /* setting the initial conditions */
    
    for (x = 0; x < 501; ++x)
      {
          p[x] = 20653.3;
          d[x] = 10;
          v[x] = 50;
       
          energy[x] = (p[x])/((d[x])*(y-1));
          e[x] = energy[x]*d[x] + 0.5*v[x]*v[x]*d[x];
      }
    for (x = 501; x < N; ++x)
      {
          p[x] = 619.6;
          d[x] = 30;
          v[x] = -50;
          
          energy[x] = (p[x])/((d[x])*(y-1));
          e[x] = energy[x]*d[x] + 0.5*v[x]*v[x]*d[x];
      }
    
    /* start of the while loop to loop through time */
    
    while (t < 0.08){
        
        /* start of the for loop to loop through the x positions */
        
        for (x = 0; x < N-1; x++)
        {
            f(d[x+1], v[x+1], e[x+1], p[x+1], fr);
            f(d[x], v[x], e[x], p[x], fl);
            
            q(d[x], v[x], e[x], ql);
            q(d[x+1], v[x+1], e[x+1], qr);
            
            cl = cs(p[x], d[x]);
            cr = cs(p[x+1], d[x+1]);
    
            pres = pressure_estimate(p[x], p[x+1], v[x], v[x+1], d[x], d[x+1], cl, cr);
            
            /* decision tree to work out the left and right wave speeds */
            
            if (pres <= p[x]){
                sl = v[x] - cl;
            }
            else {
                sl = v[x] - shock_wave(cl, pres, p[x]);
            }

            if (pres <= p[x+1]){
                sr = v[x+1] + cr;
            }
            
            else {
                sr = v[x+1] + shock_wave(cr, pres, p[x+1]);
            }
            
            /* function to work out the speed of the discontinous zone*/
            
            smid = (p[x+1] - p[x] + d[x]*v[x]*(sl - v[x]) - d[x+1]*v[x+1]*(sr - v[x+1]))/(d[x]*(sl - v[x]) - d[x+1]*(sr - v[x+1]));

            /* decision tree to work out the fluxes in between the zones. If both waves are travelling slower than the left wave the middle flux is equal to the left zone, if both waves are travelling right the middle flux is equal to the right zone. If the waves are travelling within the discontinous zone the flux is based off the wave speeds, left OR right fluxes, the speed of the discontinous zone and the fluid state in between the zones being investigated */
            
            if (sl >= 0)
            {
                f(d[x], v[x], e[x], p[x], fux);

            }

            else if (0 >= sr)
            {
                f(d[x+1], v[x+1], e[x+1], p[x+1], fux);
            }

            else if (sl <= 0 && smid >= 0)
            {
                
                qlm[0] = qmid1(d[x], v[x], sl, smid);
                qlm[1] = qmid2(d[x], v[x], sl, smid);
                qlm[2] = qmid3(d[x], v[x], smid, sl, e[x], p[x]);

                fux[0] = fl[0] + sl*(qlm[0] - ql[0]);
                fux[1] = fl[1] + sl*(qlm[1] - ql[1]);
                fux[2] = fl[2] + sl*(qlm[2] - ql[2]);
            }

            else if (smid <= 0 && sr >= 0)
            {
               
                qrm[0] = qmid1(d[x+1], v[x+1], sr, smid);
                qrm[1] = qmid2(d[x+1], v[x+1], sr, smid);
                qrm[2] = qmid3(d[x+1], v[x+1], smid, sr, e[x+1], p[x+1]);

                fux[0] = fr[0] + sr*(qrm[0] - qr[0]);
                fux[1] = fr[1] + sr*(qrm[1] - qr[1]);
                fux[2] = fr[2] + sr*(qrm[2] - qr[2]);
            }
            
            flux[0][x] = fux[0];
            flux[1][x] = fux[1];
            flux[2][x] = fux[2];
            
            qmid[0][x] = ql[0];
            qmid[1][x] = ql[1];
            qmid[2][x] = ql[2];
        }
        
        /* for loop to calculate the new fluid states from the fluxes of the half states*/
        
        for (x = 1; x < N-1; x++)
        {
            q1[x] = qhfinal(dt, flux[0][x], flux[0][x-1], qmid[0][x], dx);
            q2[x] = qhfinal(dt, flux[1][x], flux[1][x-1], qmid[1][x], dx);
            q3[x] = qhfinal(dt, flux[2][x], flux[2][x-1], qmid[2][x], dx);
        }
        
        /* for loop to calculate the new density, energy, pressure and velocity values aswell as asign them to the correcta arrays*/
        
        for (x = 1; x < N-1; ++x)
          {
              d[x] = q1[x];
              e[x] = q3[x];
              v[x] = q2[x]/q1[x];
              energy[x] = total_e(e[x], d[x], v[x]);
              p[x] = (y-1)*d[x]*energy[x];

          }
        /* calulcated the time step */
        dt = 0.9*(0.01)/(max(v, N) + sqrt(y*(max(p, N)/min(d, N))));
        t = t + dt;
    }
    for (x = 1; x < N-1; ++x)
    {
        /* writing the array to a text file */
        fprintf(data,"%d %f %f %f %f\n",x, d[x], v[x], p[x], energy[x]);
    }
    /* freeing the memory dedicated to the arrays and closing the file pointer */
    free(p);
    free(d);
    free(v);
    free(e);
    free(q1);
    free(q2);
    free(q3);
    free(q1n);
    free(q2n);
    free(q3n);
    fclose(data);
}
