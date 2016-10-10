
#include"non_uniform_grid.h"

//extern OctCell (*bodygrid)[Nx+2];
extern OctCell *bodygrid;
//extern double hx, hy;
//初始网格的生成

void load_grid();

typedef struct gridcell
{
	int pn[4];
}GRIDCELL;

//GRIDCELL (*bgrid)[Nx4grd+2];
GRIDCELL (*bgrid)[Nx];

PXYZ *gpoint;
GRIDCELL *gcell;

#define EQUAL_SPACE_STEP

void background_grid()
{
    int j,k, m;
    unsigned NO;

    double xMin, yMin;

    OctCell *pc;
    pc=bodygrid+1;
    NO=1u;
#ifndef EQUAL_SPACE_STEP
    bgrid = new GRIDCELL [Ny][Nx];
    gpoint = new PXYZ [NPOINT+1];
    gcell = new GRIDCELL [NCELL+1];	
    load_grid();
#endif        
    xMin=xLeft-2.0*hx;  //等距网格
    yMin=yLow- 2.0*hy;

    for(j=1;j<=Ny;j++) {
        for(k=1;k<=Nx;k++) {	
            //	等距网格
#ifdef  EQUAL_SPACE_STEP            				       		       
            bodygrid[(j-1)*Nx +k].xc1=hx/2.0+(k-1)*hx+xMin;
            bodygrid[(j-1)*Nx +k].yc1=hy/2.0+(j-1)*hy+yMin;	
            bodygrid[(j-1)*Nx +k].dx=hx;
            bodygrid[(j-1)*Nx +k].dy=hy;

            //------------------------------------------------------------------				
            //bodygrid[(j-1)*Nx +k]
#else 
            bodygrid[(j-1)*Nx +k].xc1 = 0.5*(gpoint[bgrid[j-1][k-1].pn[0]].x
                +gpoint[bgrid[j-1][k-1].pn[2]].x);
            bodygrid[(j-1)*Nx +k].yc1 = 0.5*(gpoint[bgrid[j-1][k-1].pn[0]].y
                +gpoint[bgrid[j-1][k-1].pn[2]].y);

            bodygrid[(j-1)*Nx +k].dx=gpoint[bgrid[j-1][k-1].pn[2]].x-gpoint[bgrid[j-1][k-1].pn[1]].x;
            bodygrid[(j-1)*Nx +k].dy=gpoint[bgrid[j-1][k-1].pn[0]].y-gpoint[bgrid[j-1][k-1].pn[1]].y; 
#endif				 
            //=================================================================				
            pc->set_invM();

            pc->flag=0;
            pc->level=0;
            pc->reflag=0;
            pc->coflag=0;
            pc->NOparent=NO;
            pc->parent=NULL;

            for(m=0;m<4;m++)  pc->children[m]=NULL;
            pc++;
            NO++;									  

            //	bodygrid[j][k].d = sqrt(bodygrid[j][k].xc1*bodygrid[j][k].xc1
            //		+bodygrid[j][k].yc1*bodygrid[j][k].yc1) - RADII;
        }			
    }
    cout<<"total background grid: "<< NO-1<<endl;	
#ifndef EQUAL_SPACE_STEP		
    delete [] gpoint;
    delete [] gcell;
    delete [] bgrid;
#endif
}



//PXYZ gpoint[NPOINT];
//struct gridcell gcell[NCELL];


void load_grid()
{
	int iy,ix;

	ifstream infile;
	infile.open(GRID_FILE);

	for(ix=1; ix < NPOINT+1; ix++)
		infile >> gpoint[ix].x >> gpoint[ix].y;

	for(ix=1; ix < NCELL+1; ix++)
		for(iy=0; iy<4; iy++) 
			infile >> gcell[ix].pn[iy];
        
	//for(iy=1; iy<=Ny4grd+2; iy++)
	for(iy=1; iy<=Ny; iy++)
	{
	//	for(ix=1; ix<=Nx4grd+2; ix++)
	for(ix=1; ix<=Nx; ix++)
		{
			//bgrid[iy-1][ix-1]=gcell[(ix-1)*(Ny4grd+2)+Ny4grd+2-(iy-1)];
			bgrid[iy-1][ix-1]=gcell[(ix-1)*Ny+Ny-(iy-1)];
		}
	}
	infile.close();
}


