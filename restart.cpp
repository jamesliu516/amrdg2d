#include"non_uniform_grid.h"

//extern OctCell (*bodygrid)[Nx+2];

//extern REAL initQ, uFree, vFree, initP,initK,initOmg;

//void jbblToShbl(JBBL *, SHBL *);
/*
void read_sol_data(string &infn)
{
    ifstream infile;
    infile.open(infn.c_str());
    string str1;

    getline(infile, str1);
    getline(infile, str1);  	
    getline(infile, str1);

    double r1;
    JBBL jb1;
    SHBL sh1;
    OctCell *lsbl;

    for(int jy=3;jy<Ny-1;jy++)
        for(int ix=3; ix< Nx-1; ix++)
        {
            lsbl= &bodygrid[jy][ix];   	
            infile>>r1>>r1>>jb1.q>>jb1.u>>jb1.v>>jb1.p>>r1>>r1;
            jbblToShbl(&jb1,&sh1);

            lsbl->dof0[0][0]=lsbl->dof[0][0]=sh1.q;
            lsbl->dof0[1][0]=lsbl->dof[1][0]=sh1.qu;  		  
            lsbl->dof0[2][0]=lsbl->dof[2][0]=sh1.qv; 		     		  
            lsbl->dof0[3][0]=lsbl->dof[3][0]=sh1.te;	  

            for(int imm=0;imm<4;++imm)
                for(int i=1; i<=nDOF-1; i++) 
                    lsbl->dof0[imm][i]=lsbl->dof[imm][i]=0.0;
        }
}	    		         
*/

void output4restart()
{
    Node *current;
    int i;
    current = HeadListAllGrid;  
    OctCell *pcell0;
    ofstream fp1;
    fp1.open("restart_data2");
    current = HeadListAllGrid;     
    while(current != NULL)
    {       
        pcell0 = current->cell;
        for(i=0;i<4;++i) fp1<<pcell0->dof[i][0]<<"  ";
        fp1<<"\n";
        current = current->next;
    }           
    fp1.close();
}

void read4restart()
{
    Node *current;
    int i,j;
    current = HeadListAllGrid;  
    OctCell *pcell0;
    ifstream fp1;
    fp1.open("IN//restart_data");
    current = HeadListAllGrid;     
    while(current != NULL)
    {       
        pcell0 = current->cell;
        for(i=0;i<4;++i){
            fp1>>pcell0->dof[i][0];
            pcell0->dof0[i][0]=pcell0->dof[i][0];
            for(j=1;j<nDOF;++j)
                pcell0->dof[i][j]=0.0;
        }

        current = current->next;
    }           
    fp1.close();
}





