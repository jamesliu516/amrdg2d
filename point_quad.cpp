#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include"vec2d.h"
const double PI = 3.14159265358979;
using namespace std;

int main()
{
    ofstream fpxx;
    double theta=20.0;
    int n0=78;// the number of equal spatial  
    int n1=78;
    int n2=78;
    int n3=78;
    Vec2D pt0(-3.0,2.0);
    Vec2D pt1(-3.0,-2.0);
    Vec2D pt2(5.0,-2.0);
    Vec2D pt3(5.0,2.0);
    /*int n0=30;// the number of equal spatial  
    int n1=100;
    int n2=10;
    int n3=100;
    Vec2D pt0(0.0,2.0);
    Vec2D pt1(0.0,0.0);
    Vec2D pt2(8.0,0.0);
    Vec2D pt3(8.0,0.8);*/

    Vec2D nml0, nml1,nml2,nml3;
    Vec2D vec0;

    nml0=(pt1-pt0).norm();
    nml1=(pt2-pt1).norm();
    nml2=(pt3-pt2).norm();
    nml3=(pt0-pt3).norm();

    double h0=fabs(pt1-pt0)/n0;
    double h1=fabs(pt2-pt1)/n1;
    double h2=fabs(pt2-pt3)/n2;
    double h3=fabs(pt0-pt3)/n3;

    double db1;
    int tn;
    tn=n0+1+n1+n2+n3-1;
    cout<<"total point  "<<tn<<endl;

    cout<<"singular point no. "<< 0<<"  "<<n0<<"  "<< n0+n1<<"  "<<n0+n1+n2<<endl;
    fpxx.open("IN//quad01_point");


    for(int i=0; i<n0+1; ++i){
        vec0=pt0+nml0*i*h0;
        fpxx<<setw(16)<<vec0[0] <<setw(16)<<vec0[1];
        fpxx<<"\n";
    }

    for(int i=1; i<n1+1; ++i){
        vec0=pt1+nml1*i*h1;
        fpxx<<setw(16)<<vec0[0]<<setw(16)<<vec0[1];
        fpxx<<"\n";
    }

    for(int i=1; i<n2+1; ++i){
        vec0=pt2+nml2*i*h2;
        fpxx<<setw(16)<<vec0[0]<<setw(16)<<vec0[1];
        fpxx<<"\n";
    }

    for(int i=1; i<n3; ++i){
        vec0=pt3+nml3*i*h3;
        fpxx<<setw(16)<<vec0[0]<<setw(16)<<vec0[1];
        fpxx<<"\n";
    }
    fpxx.close();

    return 0;
}

    




