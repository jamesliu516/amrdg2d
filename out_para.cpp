#include<sstream>
#include"non_uniform_grid.h"

using namespace std;
extern string str0;
extern double CFL;
void out_para( const string &str123)
{
    ofstream fpxx;
    string proxfile1(FILE_PROX);
		
	string  str1=str0+proxfile1;
	
	ostringstream omess1;
	
	omess1 << str1<<str123<<"_npara";
	
	fpxx.open(omess1.str().c_str());

    fpxx<<left<<setw(35)<<"CFL:"<<setw(10)<<right<<CFL<<endl;
    fpxx<<left<<setw(35)<<"1: with limter, 0:not:"
        <<setw(10)<<right<<Limiter01<<endl;
    fpxx<<left<<setw(35)<<"M (tvb para):"<<setw(10)<<right<<M4TVBM<<endl;
    fpxx<<left<<setw(35)<<"polynomial order:"<<setw(10)<<right
        <<(nDOF==6 ? 2:1)<<endl;
    fpxx<<left<<setw(35)<<"0:reflect,1:curvature modify:"<<setw(10)
        <<right<<qvlvModify<<endl;

    fpxx<<left<<setw(35)<<"grid file:"<<setw(10)<<right<<GRID_FILE<<endl;

    fpxx<<left<<setw(35)<<"total grid point:"<<setw(10)<<right<<NPOINT<<endl;
    fpxx<<left<<setw(35)<<"total grid cell:"<<setw(10)<<right<<NCELL<<endl;
    fpxx<<left<<setw(35)<<"wall point file:"<<setw(10)<<right<<WALL_POINT_FILE <<endl;

    fpxx<<left<<setw(35)<<"wall point number:"<<setw(10)<<right<<NPOINTCIRCLE <<endl;
    fpxx<<left<<setw(35)<<"nx,ny->"<<setw(10)<<right<<Nx<<","<<Ny<<endl;

    fpxx<<left<<setw(35)<<"initial hx,hy->"<<setw(10)<<right<<hx<<","<<hy<<endl;

    fpxx<<left<<setw(35)<<"adaptive AMRP:"<<setw(10)<<right<<AMRP<<endl;
    fpxx<<left<<setw(35)<<"coarese_coe:"<<setw(10)<<right<<coarse_coe<<endl;
    fpxx<<left<<setw(35)<<"refine_coe:"<<setw(10)<<right<<refine_coe<<endl;

    fpxx<<left<<setw(35)<<"coarese_coe1:"<<setw(10)<<right<<coarse_coe1<<endl;
    fpxx<<left<<setw(35)<<"refine_coe1:"<<setw(10)<<right<<refine_coe1<<endl;
    fpxx<<left<<setw(35)<<"solution adaptive number:"<<setw(10)<<right<<Nar<<endl;
    fpxx<<left<<setw(35)<<"Max adaptive number:"<<setw(10)<<right<<MaxAMR<<endl;
    fpxx<<left<<setw(35)<<"initial local refine Nr:"<<setw(10)<<right<<Nr<<endl;
    fpxx<<left<<setw(35)<<"direct local refine:"<<setw(10)<<right<<localRefine<<endl;
    fpxx<<left<<setw(35)<<"free Ma number:"<<setw(10)<<right<<FreeMa<<endl;
    
    fpxx<<left<<setw(35)<<"angle of attack:"<<setw(10)<<right<<AOA<<endl;



	fpxx.close();


}
