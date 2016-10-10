

#include"non_uniform_grid.h"
#include"vec2d.h"
#include<sstream> 

extern PXYZ Circle[NPOINTCIRCLE+1];
extern Node *HeadListAllGrid;
void formListForAllGrid();
extern PXYZ MultiCircle[N_BODY+1][500]; //points at solid wall
extern int NWallPts[N_BODY+3];
double approx_sol(double x, double y, double *doftmp, int n );
typedef struct point0123
{
    double x,y;
//	double SumVo;
//	double SumP,SumQ,SumU,SumV;//基本变量加和

//	double sxlgq,sylgq;
//	double xlgq,ylgq;

 //   JBBL jb;
}PointJBoutput;

typedef struct cell0123
{
	int cpoint[4];
}ControlCellOutput;

PointJBoutput *pointOutput=NULL;

ControlCellOutput *cellOutput=NULL;

static int NU=1,Np=1;

extern WALLP MultiCompany[][500]; //壁面输出量
extern WALLP Company[NPOINTCIRCLE+1]; //壁面输出量
void PrintNode_bk(OctCell* parent)
{
	int i,j;
	int ix,jy,kz,m;
//	double plevel;
	double x[4],y[4];//网格顶点的坐标
	double hcx,xc,yc,hcy;
	
	if(parent->children[1]==NULL)
	{   
	//	if(parent->flag>0)
		//{
		//	plevel=parent->level+1.0;		    
			//hc=h/pow(2.0,plevel);//子网格的步长	
			hcx=0.5*parent->dx;
			hcy=0.5*parent->dy;	   
			xc=parent->xc1;	    
			yc=parent->yc1;		
			//printf("xc=%15.8e,yc=%15.8e\n",parent->xc1,parent->yc1);
		//	cout<<"xc= "<<xc<<", yc= "<<yc<<endl;		
			for(m=0;m<4;m++)		
			{				
				switch(m)					
				{				
				case 0: ix=-1;jy=-1;kz=-1;break;				
				case 1: ix=1;jy=-1;kz=-1;break;				
				case 2: ix=1;jy=1;kz=-1;break;				
				case 3: ix=-1;jy=1;kz=-1;break;								
				}         
				x[m]=xc+ix*hcx;			
				y[m]=yc+jy*hcy;
		
			}
	
			for(j=0;j<4;j++)		
			{		
				if(Np==1)		
				{								
					pointOutput[Np].x=x[j];                   			
					pointOutput[Np].y=y[j];	
							
					cellOutput[NU].cpoint[j]=Np;								
					Np++; 		
				}					
				else		
				{													
					//for(i=1;i<Np;i++)	
					for(i=Np-1;i>0;i--)				
					{												
						if(fabs(x[j]-pointOutput[i].x)+fabs(y[j]-pointOutput[i].y)<1.0e-8) break;					
					//	if(i==Np-1) goto loop1;
					      if(i==1) goto loop1;				
					}								
					cellOutput[NU].cpoint[j]=i;				
					goto loop2;
loop1:
					{					
						pointOutput[Np].x=x[j];          														
						pointOutput[Np].y=y[j];					
						cellOutput[NU].cpoint[j]=Np;								
						Np++;  					
//					printf("point[%d].x=%f,point[%d].y=%f\n",Np-1,point[Np-1].x,Np-1,point[Np-1].y);		
//					getch();			
					}														    								          			
				}
loop2: 			continue;					
			}				
			NU++;		
     //   }
	}
	else 
	{
		for(int im=0;im<4;im++)
		{
			PrintNode_bk(parent->children[im]);
		}
	}
}
extern string str0;
		
void outputcell(int ij)
{
//
	ofstream fpxx;

	string proxfile1(FILE_PROX);
		
	string  str1=str0+proxfile1;
	
	ostringstream omess1;
	
	omess1 << str1 << ij<<".plt";
	
	//fpxx.open(filename1[ljm]);
	fpxx.open(omess1.str().c_str());
	
//
	long i;
	long k,j;
	
	Node *current;
	OctCell *lsbl;
	current = HeadListAllGrid;
	
	pointOutput=new PointJBoutput[NpointForOutput];
	cellOutput=new ControlCellOutput[NpointForOutput];

	while(current != NULL)
	{
		
		lsbl = current->cell;
		if(lsbl->flag % 2 == 0 && current->flg<=2 ) 
			PrintNode_bk(lsbl);
		current = current->next;
	}

//	for(i=1;i<=Nx*Ny;i++)
//	{
//		PrintNode(&bodygrid[i]);
//	}

        fpxx<<"TITLE =\"EULER SOLVER\"\n"
        <<"VARIABLES = \"X\", \"Y\"\n"
        <<"ZONE N=   "<<Np-1 << ",E=    "<<   NU-1<< ", F=FEPOINT, ET=QUADRILATERAL\n";  
        
         cout<<"Np= "<<Np-1<<", NU= "<<NU-1<<endl;
         
//	fprintf(fp,"TITLE = \" Bodygrid \"\n");
//	fprintf(fp,"VARIABLES = \"X\", \"Y\"\n");
//	fprintf(fp,"ZONE N=\t%d ,E=\t%d , F=FEPOINT, ET=QUADRILATERAL\n",Np-1,NU-1);
//	printf("Np=%d,NU=%d\n",Np,NU);

//	getch();
	for(j=1;j<Np;j++)
	{
	  fpxx<<pointOutput[j].x<<"      "<<pointOutput[j].y<<"\n";
     //   fprintf(fp,"%15.8e\t%15.8e\n",point[j].x,point[j].y);
	}
	for(k=1;k<NU;k++)
	{
		for(i=0;i<4;i++)
		{
		   //fprintf(fp,"%d\t",cell[k].cpoint[i]);
		    fpxx<<cellOutput[k].cpoint[i]<<"   ";

		}
		//fprintf(fp,"\n");
		fpxx<<"\n";
	}
	fpxx.close();
	
	delete[] pointOutput;
	delete[] cellOutput;
}

extern map<OctCell *, cellPointIndex> cellPoints;      
extern vector<pointInCell> points4out; 

void outputcell_new(int ij)//output mesh
{
//
	ofstream fpxx;

	string proxfile1(FILE_PROX);
		
	string  str1=str0+proxfile1;
	
	ostringstream omess1;
	
	omess1 << str1<<"mesh_tn" << ij<<".plt";
	
	fpxx.open(omess1.str().c_str());

        fpxx<<"TITLE =\"EULER SOLVER\"\n"
        <<"VARIABLES = \"X\", \"Y\"\n"
        <<"ZONE N=   "<<points4out.size() << ", E=    "<< cellPoints.size()<< ", F=FEPOINT, ET=QUADRILATERAL\n";  
        
	for(vector<pointInCell>::size_type j=0;j!=points4out.size();++j)
	{
	  fpxx<<points4out[j].pt[0]<<"      "<<points4out[j].pt[1]<<"\n";
	}
	for(map<OctCell *, cellPointIndex> ::iterator ik=cellPoints.begin();ik!=cellPoints.end(); ++ik)
	{
		for(int i=0;i<4;i++)
		{
		    fpxx<<ik->second.iPoint[i]<<"   ";

		}
		fpxx<<"\n";
	}
	
/*    fpxx<<"\n\n";	
    fpxx<<"TITLE =\"Wall Points\"\n"
        <<"VARIABLES = \"X\", \"Y\"\n"
        <<"ZONE I=   "<<NPOINTCIRCLE<< ", F=POINT\n";  

    for(int ii=0; ii< NPOINTCIRCLE; ++ii) {
        fpxx<<Circle[ii].x<<"  "<<Circle[ii].y<<endl;
    }	
*/	
	
	fpxx.close();
}

void computWallPresCp();
void computeCL_Cd(int ni);
//output solution
void getAllVorM();
void outputcell_sol(int ij, const string &str123)
{
	//
	ofstream fpxx;

	string proxfile1(FILE_PROX);

	string  str1=str0+proxfile1;

	ostringstream omess1;

	ostringstream omess2;
	omess2 << str1<<str123<<"cp_"<< ij<<".plt";
	omess1 << str1<<str123<<setw(6)<<setfill('0')<< ij<<".plt";

	fpxx.open(omess1.str().c_str());
double tm_tdiv;

	double u,v,Et, p, den, ma, u1,v1,p1, den1;
	double volsum, vol, entrpy,erss, tp, tper;
    double rtp000;
    
    if(isOutVorM) getAllVorM();

	fpxx<<"TITLE =\"EULER SOLVER\"\n"
		<< "VARIABLES = \"X\", \"Y\",   \"DEN\" ,\"U\",\"V\", \"PRESS\" ,\"MA\", \"xingqu\"\n"
		<<"ZONE N=   "<<points4out.size() << ", E=    "<< cellPoints.size()<< ", F=FEPOINT, ET=QUADRILATERAL\n";  

	for(vector<pointInCell>::size_type j=0;j!=points4out.size();++j)
	{
		volsum=0.0;
        tm_tdiv=0.0;
		u=v=p=den=0.0;
		//  if(j==0){int iitm; cout<<"look look"<< points4out[j].nmc<<endl; cin>>iitm;}

		for(int ii=0; ii< points4out[j].nmc; ++ii){
			vol=1.0/((points4out[j].mcell[ii])->dx*(points4out[j].mcell[ii])->dy);
			volsum+=vol;
			den1=(points4out[j].mcell[ii])->dof[0][0];	     	     
			u1=(points4out[j].mcell[ii])->dof[1][0]/den1;
			v1=(points4out[j].mcell[ii])->dof[2][0]/den1;
			Et=(points4out[j].mcell[ii])->dof[3][0];	     
			p1=GAM11*(Et-0.5*den1*(u1*u1+v1*v1));

			den+=vol*den1;
			u+=vol*u1;
			v+=vol*v1;
			p+=vol*p1;
           // rtp000=(((points4out[j].mcell[ii])->trb==true)?1.0:0.0);
            rtp000=(points4out[j].mcell[ii])->kxrcf;
         //   rtp000=(points4out[j].mcell[ii])->vorM;
       //     tm_tdiv+=vol*(points4out[j].mcell[ii])->tdiv;
            tm_tdiv+=vol*rtp000;
		}

		u/=volsum;
		v/=volsum;
		den/=volsum;
		p/=volsum;
        tm_tdiv/=volsum;

		entrpy = p/pow(den,GAMMA);
		erss=fabs(initP/pow(initQ,GAMMA)-entrpy)/(initP/pow(initQ,GAMMA));
		tp = 0.5 * den * (u *u + v * v) + p;
		tper = fabs(tp-0.5*initQ*(uFree*uFree + vFree*vFree)-initP);

		ma=sqrt(u*u+v*v)/sqrt(GAMMA*p/den);

		fpxx<<fixed<<points4out[j].pt[0]<<"  "<<points4out[j].pt[1]<<"  "<<setw(15)<<den<<" " << setw(15)<<u <<" " << setw(15)<<v <<" "
			<< setw(15)<< p << " "<<setw(15)<<ma<<" "
			<<"  "<<setw(15)<<tm_tdiv << '\n'; 
	}
	for(map<OctCell *, cellPointIndex> ::iterator ik=cellPoints.begin();ik!=cellPoints.end(); ++ik)
	{
		for(int i=0;i<4;i++)
		{
			fpxx<<ik->second.iPoint[i]<<"   ";

		}
		fpxx<<"\n";
	}

    for (int jj=0; jj<N_BODY;++jj){
        fpxx<<"\n\n";	
        fpxx<<"TITLE =\"Wall Points\"\n"
            // <<"VARIABLES = \"X\", \"Y\"\n"
            << "VARIABLES = \"X\", \"Y\",   \"DEN\" ,\"U\",\"V\", \"PRESS\" ,\"MA\",\"xingqu\"\n"
            <<"ZONE I=   "<<NWallPts[jj]<< ", F=POINT\n";  

        for(int ii=0; ii< NWallPts[jj]; ++ii) {
            fpxx<<MultiCircle[jj][ii].x<<"  "<<MultiCircle[jj][ii].y<<"  "<<0.0<<"  "<<0.0
                <<"  "<<0.0<<"  "<<0.0<<"  "<<0.0<<"  "<<0.0<<"\n";
        }	
    }

    computWallPresCp();
    fpxx.close();
    fpxx.clear();
    fpxx.open(omess2.str().c_str());
    for (int jj=0; jj<N_BODY;++jj){
        fpxx<<"TITLE =\"Cp\"\n"
            // <<"VARIABLES = \"X\", \"Y\"\n"
            << "VARIABLES = \"X\",\"Cp\"\n"
            <<"ZONE I=   "<<NWallPts[jj]<< ", F=POINT\n";  

        for(int ii=0; ii< NWallPts[jj]; ++ii) {
            fpxx<<MultiCircle[jj][ii].x<<"  "<<MultiCompany[jj][ii].cp<<"\n";
        }	
    }
    fpxx.close();	
    fpxx.clear();	
    computeCL_Cd(ij);
}

void output4restart()
{
    Node *current;
    int i;
    current = HeadListAllGrid;  
    OctCell *pcell0;
    ofstream fp1;

	string  str1=str0+"restart_data";
    fp1.open(str1.c_str());
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





