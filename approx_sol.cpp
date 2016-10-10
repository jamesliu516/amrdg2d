
#include"vec2d.h"

//x,y belong to [-1,1]X[-1,1] 需要将  [x_{i-1/2} x_{i+1/2}] => [-1,1]
double approx_sol(double x, double y, double *doftmp, int n )
{

    double uh;
    if(n==6) {
      uh=doftmp[0]+doftmp[1]*x+doftmp[2]*y+doftmp[3]*x*y
       +doftmp[4]*(x*x-1.0/3.0) + doftmp[5]*(y*y-1.0/3.0); 
     }
     else if(n==3) {
      uh=doftmp[0]+doftmp[1]*x+doftmp[2]*y;
     }
     else if(n==1) {
        uh=doftmp[0];
     }
       
    return uh;
}

//直接x,y是原始计算区域的坐标 [x_{i-1/2} x_{i+1/2}] 
double n_approx_sol(double x,double y, double xc1, double yc1, double dx, double dy, double *doftmp, int n )
{
    double uh;
    if(n==6) {
      uh=doftmp[0]+doftmp[1]*2.0*(x-xc1)/dx+doftmp[2]*2.0*(y-yc1)/dy+doftmp[3]*4.0*(x-xc1)*(y-yc1)/(dx*dy)
       +doftmp[4]*((x-xc1)*(x-xc1)*4.0/(dx*dx)-1.0/3.0) + doftmp[5]*((y-yc1)*(y-yc1)*4.0/(dy*dy)-1.0/3.0); 
     }
     else if(n==3) {
      uh=doftmp[0]+doftmp[1]*2.0*(x-xc1)/dx+doftmp[2]*2.0*(y-yc1)/dy;
     }
     else if(n==1) {
        uh=doftmp[0];
     }
       
    return uh;
}



double n_approx_sol(const double xi[2],  const double xc[2], double dx, double dy, double *doftmp, int n )
{
    double uh;
    if(n==6) {
      uh=doftmp[0]+doftmp[1]*2.0*(xi[0]-xc[0])/dx+doftmp[2]*2.0*(xi[1]-xc[1])/dy+doftmp[3]*4.0*(xi[0]-xc[0])*(xi[1]-xc[1])/(dx*dy)
       +doftmp[4]*((xi[0]-xc[0])*(xi[0]-xc[0])*4.0/(dx*dx)-1.0/3.0) + doftmp[5]*((xi[1]-xc[1])*(xi[1]-xc[1])*4.0/(dy*dy)-1.0/3.0); 
     }
     else if(n==3) {
      uh=doftmp[0]+doftmp[1]*2.0*(xi[0]-xc[0])/dx+doftmp[2]*2.0*(xi[1]-xc[1])/dy;
     }
     else if(n==1) {
        uh=doftmp[0];
     }       
    return uh;
}

