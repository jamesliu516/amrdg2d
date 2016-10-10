#include"non_uniform_grid.h"

int newFlagpoint(double x, double y);////应用射线方法
void findrefp(OctCell *pp[], OctCell *parent);

extern OctCell *bodygrid;

//extern OctCell *refcell[][NneedRefineGrid];
extern vector<OctCell *> refcell[Nr+1];

static int re=0;  //for creat grid at beginning
//int kref=0;//for creat grid at beginning

//extern re, kref;

extern PXYZ Circle[NPOINTCIRCLE+1];

extern PXYZ MultiCircle[N_BODY+1][500]; //points at solid wall
extern int NWallPts[N_BODY+3];
extern PXYZ ext_wall[NPT_EXT_WALL+1]; //points at external computational domain
int newFlagpoint(double x, double y, int ,PXYZ waoint[] );
//网格单元类型的判定:0流体网格单元, -1,-物体内部非计算网格单元, 2-与物体相交的网格单元;
void charggrid(OctCell *pg)//pg为网格单元的指针
{
    int flag[N_BODY+1],i;
    int exflag;
    //double plevel;
    double x[4],y[4];//网格顶点的坐标
    double hcx,xc,yc, hcy;
    int pflag[N_BODY+1][4];
    int n0,nm1,n2; 
    //plevel=pg->level+1.0;	
    //hc=h/pow(2.0,plevel);//子网格的步长	
    hcx=0.5*pg->dx;
    hcy=0.5*pg->dy;

    xc=pg->xc1;
    yc=pg->yc1;
#if IS_TUBE ==1
    exflag=newFlagpoint(xc,yc,NPT_EXT_WALL,ext_wall);
    pg->exflag=(exflag==0?-1:0);
#endif

    x[0]=xc-hcx; y[0]=yc-hcy;
    x[1]=xc-hcx; y[1]=yc+hcy;	
    x[2]=xc+hcx; y[2]=yc-hcy;
    x[3]=xc+hcx; y[3]=yc+hcy;
    //-------------------
    for (int jj=0; jj<N_BODY;++jj){
        for(i=0;i<4;i++) {
            pflag[jj][i]=newFlagpoint(x[i],y[i],
                NWallPts[jj], MultiCircle[jj]);
        }

        if(pflag[jj][0]==-1 && pflag[jj][1]==-1
            && pflag[jj][2]==-1 && pflag[jj][3]==-1)
        {
            flag[jj]=-1;
            // pg->flag=-1;
        }
        else if(pflag[jj][0]==0 && pflag[jj][1]==0
            && pflag[jj][2]==0 && pflag[jj][3]==0)
        {
            flag[jj]=0;
            //  pg->flag=0;
        }
        else 
        {			
            flag[jj]=2;			
            //     pg->flag=2;
            //  	printf("xc=%f,yc=%f\n",xc,yc);
        }
    }

    n0=nm1=n2=0;
    
    for(int jj=0;jj<N_BODY;++jj){
        if(flag[jj]==0) ++n0;
        if(flag[jj]==2) ++n2;
        if(flag[jj]==-1) ++nm1;
    }

    if(n0==N_BODY) pg->flag=0;
    else if(n2>0) pg->flag=2;
    else pg->flag=-1;

    //如果该网格单元需要加密则将其地址存入数组中
    if(pg->flag==2)
    {
        //refcell[re][kref]=pg;
        //kref++;
        refcell[re].push_back(pg);
    }
}
//only check the flag of the cell, do not add the cell to refcell
void charggrid0(OctCell *pg)//pg为网格单元的指针
{
    int flag[N_BODY+1],i;
    int exflag;
    //double plevel;
    double x[4],y[4];//网格顶点的坐标
    double hcx,xc,yc, hcy;
    int pflag[N_BODY+1][4];
    int n0,nm1,n2; 
    //plevel=pg->level+1.0;	
    //hc=h/pow(2.0,plevel);//子网格的步长	
    hcx=0.5*pg->dx;
    hcy=0.5*pg->dy;

    xc=pg->xc1;
    yc=pg->yc1;
#if IS_TUBE ==1
    exflag=newFlagpoint(xc,yc,NPT_EXT_WALL,ext_wall);
    pg->exflag=(exflag==0?-1:0);
#endif

    x[0]=xc-hcx; y[0]=yc-hcy;
    x[1]=xc-hcx; y[1]=yc+hcy;	
    x[2]=xc+hcx; y[2]=yc-hcy;
    x[3]=xc+hcx; y[3]=yc+hcy;
    //-------------------
    for (int jj=0; jj<N_BODY;++jj){
        for(i=0;i<4;i++) {
            pflag[jj][i]=newFlagpoint(x[i],y[i],
                NWallPts[jj], MultiCircle[jj]);
        }

        if(pflag[jj][0]==-1 && pflag[jj][1]==-1
            && pflag[jj][2]==-1 && pflag[jj][3]==-1)
        {
            flag[jj]=-1;
            // pg->flag=-1;
        }
        else if(pflag[jj][0]==0 && pflag[jj][1]==0
            && pflag[jj][2]==0 && pflag[jj][3]==0)
        {
            flag[jj]=0;
            //  pg->flag=0;
        }
        else 
        {			
            flag[jj]=2;			
            //     pg->flag=2;
            //  	printf("xc=%f,yc=%f\n",xc,yc);
        }
    }

    n0=nm1=n2=0;
    
    for(int jj=0;jj<N_BODY;++jj){
        if(flag[jj]==0) ++n0;
        if(flag[jj]==2) ++n2;
        if(flag[jj]==-1) ++nm1;
    }

    if(n0==N_BODY) pg->flag=0;
    else if(n2>0) pg->flag=2;
    else pg->flag=-1;

}
/*
void charggrid0(OctCell *pg)//pg为网格单元的指针
{
	int flag,i,exflag;
	//double plevel;
	double x[4],y[4];//网格顶点的坐标
	double hcx,xc,yc, hcy;
       int pflag[4];
	//plevel=pg->level+1.0;	
	//hc=h/pow(2.0,plevel);//子网格的步长	
	hcx=0.5*pg->dx;
	hcy=0.5*pg->dy;

	xc=pg->xc1;
	yc=pg->yc1;
    exflag=newFlagpoint(xc,yc,NPT_EXT_WALL,ext_wall);
    pg->exflag=(exflag==0?-1:0);
	x[0]=xc-hcx; y[0]=yc-hcy;
	x[1]=xc-hcx; y[1]=yc+hcy;	
	x[2]=xc+hcx; y[2]=yc-hcy;
	x[3]=xc+hcx; y[3]=yc+hcy;

	for(i=0;i<4;i++) {
          pflag[i]=newFlagpoint(x[i],y[i]);
	}
	
       if(pflag[0]==-1 && pflag[1]==-1 && pflag[2]==-1 && pflag[3]==-1)
	{
		flag=-1;
		pg->flag=-1;
	}
	else if(pflag[0]==0 && pflag[1]==0 && pflag[2]==0 && pflag[3]==0)
	{
		flag=0;
		pg->flag=0;
	}
	else 
	{			
		flag=2;			
		pg->flag=2;
  //  	printf("xc=%f,yc=%f\n",xc,yc);
	}
    //如果该网格单元需要加密则将其地址存入数组中
}
*/
/*---对父亲节点所对应的孩子节点赋值---*/
void childrenvalue(OctCell *p1, OctCell *parent,int m,double hcx, double hcy)
{      
	int i,j;
	switch(m)
	{	
	case 0: i=-1;j=-1;break;
	case 1: i=1;j=-1;break;
	case 2: i=1;j=1;break;
	case 3: i=-1;j=1;break;	
	}	   
//	   p1->flag=0;
	   p1->level=parent->level+1;
          p1->NOparent=parent->NOparent;
	   p1->NOchildren=m;
          p1->level0=parent->level0;
	//   p1->reflag=0; //constructor can do this
	//   p1->coflag=0;
	   p1->parent=parent;
	   p1->dx=hcx;
	   p1->dy=hcy;
	   p1->set_invM();
	   
//	   for(int n=0;n<4;n++)	
//	   {                    
//		   p1->children[n]=NULL;			
//	   }	
	   							
	   p1->xc1=parent->xc1+i*hcx*0.5;				
	   p1->yc1=parent->yc1+j*hcy*0.5;	
	   			
}

void childrenvalue(OctCell *p1, OctCell *parent,
    double hcx, double hcy, int m)
{      
	int i,j;
	switch(m)
	{	
	case 0: i=-1;j=-1;break;
	case 1: i=1;j=-1;break;
	case 2: i=1;j=1;break;
	case 3: i=-1;j=1;break;	
	}	   
//	   p1->flag=0;
	   p1->level=parent->level+1;
          p1->NOparent=parent->NOparent;
          p1->level0=parent->level0;
	   p1->NOchildren=m;
	//   p1->reflag=0; //constructor can do this
	//   p1->coflag=0;
	   p1->parent=parent;
	   p1->dx=hcx;
	   p1->dy=hcy;
	   p1->set_invM();
	   
//	   for(int n=0;n<4;n++)	
//	   {                    
//		   p1->children[n]=NULL;			
//	   }	
	   							
	   p1->xc1=parent->xc1+i*hcx*0.5;				
	   p1->yc1=parent->yc1+j*hcy*0.5;	
	   			
}

//添加子网格单元，并对子网格单元赋值
void creatchildren( OctCell *parent)
{
	int m,preflag;
	//double plevel; 
	double hcx, hcy;//子网格的步长
	OctCell *children[4];

	preflag=parent->reflag;

	if(preflag==1) return;//判断该网格单元是否加密
	else 
	{
	//	plevel=parent->level+1;	

		//hc=h/pow(2.0,plevel);//子网格的步长	
		hcx=0.5*parent->dx;
		hcy=0.5*parent->dy;

		for(m=0;m<4;m++)
		{		
			//children[m]=(struct bodygrid *)malloc(LEN);
			children[m]=new OctCell(hcx,hcy);	
			parent->children[m]=children[m];
			childrenvalue(children[m],parent,m,hcx, hcy);			
			charggrid(children[m]);	
		} 
	}
}

void creatgrid()
{   
   int ir;	
   int ic,Nc;
   int kr,Nk;
   
   int i;
   
   OctCell *parent;
   
   OctCell *pp[48]={NULL};//物面相交的网格邻近的需要加密的网格单元地址  
   Nc=Nx*Ny;
 //  kref=0; //
   re=0;
   for(ir=0;ir<Nr;ir++)
   {
	   if(ir==0)
	   {
		   for(ic=1;ic<=Nc;ic++)//找出z=0截面上需要加密的网格单元
		   {			               
			   charggrid(&bodygrid[ic]);
		   }
	   }
	   
	   Nk=refcell[re].size()-1;
	//   Nk=kref-1;
	//   kref=0;
	   re++;
	   
	   for(kr=0;kr<=Nk;kr++)//Nk与物面相交的网格数(ljm Nk+1)
	   {
	//	   OctCell *pp[48];//物面相交的网格邻近的需要加密的网格单元地址  
               //  parent=refcell[re-1][kr];
                 parent=refcell[re-1][kr];
		   findrefp(pp,parent);	//找到物面相交的网格邻近的需要加密的网格单元           
		   creatchildren(parent);
		   parent->reflag=1;
		   for(i=0;i<48;i++)
		   {                
		      creatchildren(pp[i]);
		       pp[i]->reflag=1;
		   }
	   }
   }
   
   for(int itm=0; itm<Nr+1; ++itm) 
        refcell[itm].clear();
}



//parent为叶子节点
void refineOneCell(OctCell *parent, int rfN)
{
	double hcx,hcy;//子网格的步长
	OctCell *children[4]={NULL};
//	double plevel;
	int m;

//	plevel=parent->level+1;	
//	hc=h/pow(2.0,plevel);//子网格的步长	

    if( rfN >= 1)
	{
		parent->reflag=1; /*1月2号修改*/
		hcx=0.5*parent->dx;
		hcy=0.5*parent->dy;
		for(m=0;m<4;m++)
		{		
			children[m]=new OctCell(hcx,hcy);	
			parent->children[m]=children[m];
			childrenvalue(children[m],parent,m,hcx, hcy);			
			refineOneCell(children[m], rfN-1);
		}
	}
	else return;
}

//指定区域的加密
void refineCicleDemain(int rfN)// 1,2,3,4,5次加密
{
	int i;
	int jk1,jk2,jk3,jk4,jk5;

	switch(rfN)
	{	
	case 5: 
		for(i=1;i<=Nx*Ny;i++) 
		{			
			if(jk1=(bodygrid[i].xc1 >= Ellipse_x - 2 * hx && bodygrid[i].xc1 <= Ellipse_x1 + 2 * hx		    
				&& bodygrid[i].yc1 >= Ellipse_y - 2 * hy && bodygrid[i].yc1 <= Ellipse_y1 + 2 * hy))		
			{					
				refineOneCell(&bodygrid[i], 5);
			}		
			else if(jk2=((!jk1) && bodygrid[i].xc1 > Ellipse_x - 4*hx && bodygrid[i].xc1 <= Ellipse_x1 + 4*hx 					
				&& bodygrid[i].yc1 >= Ellipse_y-4*hy && bodygrid[i].yc1 <= Ellipse_y1+4*hy))	
			{					
				refineOneCell(&bodygrid[i], 4);						
			}	    
			else if(jk3=((!jk2) && bodygrid[i].xc1 > Ellipse_x-6*hx && bodygrid[i].xc1 <= Ellipse_x1 + 6*hx 		   
				&& bodygrid[i].yc1 >=  Ellipse_y-6*hy && bodygrid[i].yc1 < Ellipse_y1 + 6*hy ) )
			{					
				refineOneCell(&bodygrid[i], 3);		
			}
			else if(jk4=((!jk3) && bodygrid[i].xc1 > Ellipse_x-8*hx && bodygrid[i].xc1 <= Ellipse_x1 + 8*hx 		   
				&& bodygrid[i].yc1 >=  Ellipse_y-8*hy && bodygrid[i].yc1 < Ellipse_y1 + 8*hy) )
			{					
				refineOneCell(&bodygrid[i], 2);		
			}				
			else if(jk5=((!jk4) && bodygrid[i].xc1 > Ellipse_x-10*hx && bodygrid[i].xc1 <= Ellipse_x1 + 10*hx 		   
				&& bodygrid[i].yc1 >=  Ellipse_y-10*hy && bodygrid[i].yc1 < Ellipse_y1 + 10*hy) )
			{					
				refineOneCell(&bodygrid[i], 1);		
			}				
		}
		break;

	case 4:
		for(i=1;i<=Nx*Ny;i++) 
		{			
			if(jk1=(bodygrid[i].xc1 >= Ellipse_x - 2 * hx && bodygrid[i].xc1 <= Ellipse_x1+ 2 * hx		    
				&& bodygrid[i].yc1 >= Ellipse_y - 2 * hy && bodygrid[i].yc1 <= Ellipse_y1 + 2 * hy))		
			{					
				refineOneCell(&bodygrid[i], 4);
			}		
			else if(jk2=((!jk1) && bodygrid[i].xc1 > Ellipse_x - 4*hx && bodygrid[i].xc1 <= Ellipse_x1 + 4*hx 					
				&& bodygrid[i].yc1 >= Ellipse_y-4*hy && bodygrid[i].yc1 <= Ellipse_y1+4*hy))	
			{					
				refineOneCell(&bodygrid[i], 3);						
			}	    
			else if(jk3=((!jk2) && bodygrid[i].xc1 > Ellipse_x-6*hx && bodygrid[i].xc1 <= Ellipse_x1+ 6*hx 		   
				&& bodygrid[i].yc1 >=  Ellipse_y-6*hy && bodygrid[i].yc1 < Ellipse_y1 + 6*hy ) )
			{					
				refineOneCell(&bodygrid[i], 2);		
			}
			else if(jk4=((!jk3) && bodygrid[i].xc1 > Ellipse_x-8*hx && bodygrid[i].xc1 <= Ellipse_x1 + 8*hx 		   
				&& bodygrid[i].yc1 >=  Ellipse_y-8*hy && bodygrid[i].yc1 < Ellipse_y1 + 8*hy) )
			{					
				refineOneCell(&bodygrid[i], 1);		
			}							
		}
		break;
	case 3:
		for(i=1;i<=Nx*Ny;i++) 
		{			
			if(jk1=(bodygrid[i].xc1 >= Ellipse_x - 2 * hx && bodygrid[i].xc1 <= Ellipse_x1 + 2 * hx		    
				&& bodygrid[i].yc1 >= Ellipse_y - 2 * hy && bodygrid[i].yc1 <= Ellipse_y1 + 2 * hy))		
			{					
				refineOneCell(&bodygrid[i], 3);
			}		
			else if(jk2=((!jk1) && bodygrid[i].xc1 > Ellipse_x - 4*hx && bodygrid[i].xc1 <= Ellipse_x1 + 4*hx 					
				&& bodygrid[i].yc1 >= Ellipse_y-4*hy && bodygrid[i].yc1 <= Ellipse_y1+4*hy))	
			{					
				refineOneCell(&bodygrid[i], 2);						
			}	    
			else if(jk3=((!jk2) && bodygrid[i].xc1 > Ellipse_x-6*hx && bodygrid[i].xc1 <= Ellipse_x1 + 6*hx 		   
				&& bodygrid[i].yc1 >=  Ellipse_y-6*hy && bodygrid[i].yc1 < Ellipse_y1 + 6*hy ) )
			{					
				refineOneCell(&bodygrid[i], 1);		
			}						
		}
		break;
	case 2:
		for(i=1;i<=Nx*Ny;i++) 
		{			
			if(jk1=(bodygrid[i].xc1 >= Ellipse_x - 2 * hx && bodygrid[i].xc1 <= Ellipse_x1 + 2 * hx		    
				&& bodygrid[i].yc1 >= Ellipse_y - 2 * hy && bodygrid[i].yc1 <= Ellipse_y1+ 2 * hy))		
			{					
				refineOneCell(&bodygrid[i], 2);
			}		
			else if(jk2=((!jk1) && bodygrid[i].xc1 > Ellipse_x - 4*hx && bodygrid[i].xc1 <= Ellipse_x1+ 4*hx 					
				&& bodygrid[i].yc1 >= Ellipse_y-4*hy && bodygrid[i].yc1 <= Ellipse_y1+4*hy))	
			{					
				refineOneCell(&bodygrid[i], 1);						
			}	    
		}
		break;
	case 1:
		for(i=1;i<=Nx*Ny;i++) 
		{			
			if(jk1=(bodygrid[i].xc1 >= Ellipse_x - 2 * hx && bodygrid[i].xc1 <= Ellipse_x1+ 2 * hx		    
				&& bodygrid[i].yc1 >= Ellipse_y - 2 * hy && bodygrid[i].yc1 <= Ellipse_y1+ 2 * hy))		
			{					
				refineOneCell(&bodygrid[i], 1);
			}		
		}
		break;
	default:
		cout<<"Only refine 1,2,3,4,5 time.\n";
		cout<<"You input: "<<rfN<<". No refine out domain.\n";
		break;
	}
}
/*
//指定区域的加密
void refineCicleDemain(int rfN)// 1,2,3,4,5次加密
{
	int i;
	int jk1,jk2,jk3,jk4,jk5;

	switch(rfN)
	{	
	case 5: 
		for(i=1;i<=Nx*Ny;i++) 
		{			
			if(jk1=(bodygrid[i].xc1 >= -Ellipse_x - 2 * hx && bodygrid[i].xc1 <= Ellipse_x + 2 * hx		    
				&& bodygrid[i].yc1 >= -Ellipse_y - 2 * hy && bodygrid[i].yc1 <= Ellipse_y + 2 * hy))		
			{					
				refineOneCell(&bodygrid[i], 5);
			}		
			else if(jk2=((!jk1) && bodygrid[i].xc1 > -Ellipse_x - 4*hx && bodygrid[i].xc1 <= Ellipse_x + 4*hx 					
				&& bodygrid[i].yc1 >= -Ellipse_y-4*hy && bodygrid[i].yc1 <= Ellipse_y+4*hy))	
			{					
				refineOneCell(&bodygrid[i], 4);						
			}	    
			else if(jk3=((!jk2) && bodygrid[i].xc1 > -Ellipse_x-6*hx && bodygrid[i].xc1 <= Ellipse_x + 6*hx 		   
				&& bodygrid[i].yc1 >=  -Ellipse_y-6*hy && bodygrid[i].yc1 < Ellipse_y + 6*hy ) )
			{					
				refineOneCell(&bodygrid[i], 3);		
			}
			else if(jk4=((!jk3) && bodygrid[i].xc1 > -Ellipse_x-8*hx && bodygrid[i].xc1 <= Ellipse_x + 8*hx 		   
				&& bodygrid[i].yc1 >=  -Ellipse_y-8*hy && bodygrid[i].yc1 < Ellipse_y + 8*hy) )
			{					
				refineOneCell(&bodygrid[i], 2);		
			}				
			else if(jk5=((!jk4) && bodygrid[i].xc1 > -Ellipse_x-10*hx && bodygrid[i].xc1 <= Ellipse_x + 10*hx 		   
				&& bodygrid[i].yc1 >=  -Ellipse_y-10*hy && bodygrid[i].yc1 < Ellipse_y + 10*hy) )
			{					
				refineOneCell(&bodygrid[i], 1);		
			}				
		}
		break;

	case 4:
		for(i=1;i<=Nx*Ny;i++) 
		{			
			if(jk1=(bodygrid[i].xc1 >= -Ellipse_x - 2 * hx && bodygrid[i].xc1 <= Ellipse_x+ 2 * hx		    
				&& bodygrid[i].yc1 >= -Ellipse_y - 2 * hy && bodygrid[i].yc1 <= Ellipse_y + 2 * hy))		
			{					
				refineOneCell(&bodygrid[i], 4);
			}		
			else if(jk2=((!jk1) && bodygrid[i].xc1 > -Ellipse_x - 4*hx && bodygrid[i].xc1 <= Ellipse_x + 4*hx 					
				&& bodygrid[i].yc1 >= -Ellipse_y-4*hy && bodygrid[i].yc1 <= Ellipse_y+4*hy))	
			{					
				refineOneCell(&bodygrid[i], 3);						
			}	    
			else if(jk3=((!jk2) && bodygrid[i].xc1 > -Ellipse_x-6*hx && bodygrid[i].xc1 <= Ellipse_x+ 6*hx 		   
				&& bodygrid[i].yc1 >=  -Ellipse_y-6*hy && bodygrid[i].yc1 < Ellipse_y + 6*hy ) )
			{					
				refineOneCell(&bodygrid[i], 2);		
			}
			else if(jk4=((!jk3) && bodygrid[i].xc1 > -Ellipse_x-8*hx && bodygrid[i].xc1 <= Ellipse_x + 8*hx 		   
				&& bodygrid[i].yc1 >=  -Ellipse_y-8*hy && bodygrid[i].yc1 < Ellipse_y + 8*hy) )
			{					
				refineOneCell(&bodygrid[i], 1);		
			}							
		}
		break;
	case 3:
		for(i=1;i<=Nx*Ny;i++) 
		{			
			if(jk1=(bodygrid[i].xc1 >= -Ellipse_x - 2 * hx && bodygrid[i].xc1 <= Ellipse_x + 2 * hx		    
				&& bodygrid[i].yc1 >= -Ellipse_y - 2 * hy && bodygrid[i].yc1 <= Ellipse_y + 2 * hy))		
			{					
				refineOneCell(&bodygrid[i], 3);
			}		
			else if(jk2=((!jk1) && bodygrid[i].xc1 > -Ellipse_x - 4*hx && bodygrid[i].xc1 <= Ellipse_x + 4*hx 					
				&& bodygrid[i].yc1 >= -Ellipse_y-4*hy && bodygrid[i].yc1 <= Ellipse_y+4*hy))	
			{					
				refineOneCell(&bodygrid[i], 2);						
			}	    
			else if(jk3=((!jk2) && bodygrid[i].xc1 > -Ellipse_x-6*hx && bodygrid[i].xc1 <= Ellipse_x + 6*hx 		   
				&& bodygrid[i].yc1 >=  -Ellipse_y-6*hy && bodygrid[i].yc1 < Ellipse_y + 6*hy ) )
			{					
				refineOneCell(&bodygrid[i], 1);		
			}						
		}
		break;
	case 2:
		for(i=1;i<=Nx*Ny;i++) 
		{			
			if(jk1=(bodygrid[i].xc1 >= -Ellipse_x - 2 * hx && bodygrid[i].xc1 <= Ellipse_x + 2 * hx		    
				&& bodygrid[i].yc1 >= -Ellipse_y - 2 * hy && bodygrid[i].yc1 <= Ellipse_y + 2 * hy))		
			{					
				refineOneCell(&bodygrid[i], 2);
			}		
			else if(jk2=((!jk1) && bodygrid[i].xc1 > -Ellipse_x - 4*hx && bodygrid[i].xc1 <= Ellipse_x + 4*hx 					
				&& bodygrid[i].yc1 >= -Ellipse_y-4*hy && bodygrid[i].yc1 <= Ellipse_y+4*hy))	
			{					
				refineOneCell(&bodygrid[i], 1);						
			}	    
		}
		break;
	case 1:
		for(i=1;i<=Nx*Ny;i++) 
		{			
			if(jk1=(bodygrid[i].xc1 >= -Ellipse_x - 2 * hx && bodygrid[i].xc1 <= Ellipse_x + 2 * hx		    
				&& bodygrid[i].yc1 >= -Ellipse_y - 2 * hy && bodygrid[i].yc1 <= Ellipse_y + 2 * hy))		
			{					
				refineOneCell(&bodygrid[i], 1);
			}		
		}
		break;
	default:
		cout<<"Only refine 1,2,3,4,5 time.\n";
		cout<<"You input: "<<rfN<<". No refine out domain.\n";
		break;
	}
}

*/



