#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include"vec2d.h"
const double PI = 3.14159265358979;
using namespace std;
const double rrr=0.5;
int main()
{
    ofstream fpxx;
    double theta=20.0;
    double stp;
    int NWallPts=100;
    int n0=60;// the number of equal spatial  
    double x0,y0;
    x0=0.0;
    y0=0.0;
    Vec2D vec0;
    Vec2D pt0;
    Vec2D pt1;
    Vec2D nml0;

    int tn;
    tn=n0+NWallPts;
    cout<<"total point  "<<tn<<endl;

    cout<<"singular point no. "<< 0<<"  "<<NWallPts+1<<endl;
    fpxx.open("IN//halfCircle");

    stp = PI/NWallPts;
    for(int i=0; i <= NWallPts; i++) {
        vec0[0]=x0+rrr*cos(PI*0.5+i*stp);
        vec0[1]=y0+rrr*sin(PI*0.5+i*stp);
        fpxx<<setw(16)<<vec0[0] <<setw(16)<<vec0[1];
    }
    pt0[0]=x0+rrr*cos(1.5*PI);
    pt0[1]=y0+rrr*sin(PI*1.5);

    pt1[0]=x0+rrr*cos(0.5*PI);
    pt1[1]=y0+rrr*sin(PI*0.5);

    nml0=(pt1-pt0).norm();
    double h0=fabs(pt1-pt0)/n0;
    for(int i=1; i<n0; ++i){
        vec0=pt0+nml0*i*h0;
        fpxx<<setw(16)<<vec0[0] <<setw(16)<<vec0[1];
        fpxx<<"\n";
    }
//-----------------------------
    fpxx.close();

    return 0;
}

    




