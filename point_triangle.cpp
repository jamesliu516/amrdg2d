#include<iostream>
#include<fstream>
#include<iomanip>
#include<cmath>
const double PI = 3.14159265358979;
using namespace std;

int main()
{
    ofstream fpxx;
    double h=1.0;
    double theta=30.0;
    int nx=100;
    int ny=100;
    double hx,hy;
    double db_tmp;
    int tn;
    tn=nx+1+ny+nx-1;
    cout<<"total point  "<<tn<<endl;

    cout<<"singular point no. "<< 0<<"  "<<nx<<"  "<< nx+ny<<endl;
    fpxx.open("IN//tri_Pro");

    hx=h/nx;
    db_tmp=tan(PI*theta/180);
    hy=2.0*h*tan(PI*theta/180)/ny;

    for(int i=0; i<nx+1; ++i){
        fpxx<<setw(16)<<i*hx<<setw(16)<<-tan(PI*theta/180)*i*hx;
        fpxx<<"\n";
    }

    for(int i=0; i<ny; ++i){
        fpxx<<setw(16)<<h<<setw(16)<<-db_tmp*nx*hx+(i+1)*hy;
        fpxx<<"\n";
    }

    for(int i=nx-1;i>0;--i){
        fpxx<<setw(16)<<i*hx<<setw(16)<<db_tmp*i*hx<<"\n";
    }

    fpxx.close();

    return 0;
}

    




