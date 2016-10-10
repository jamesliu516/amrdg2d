// solution AMR
#include"non_uniform_grid.h"
#include"findneighbor.h"
extern vector<OctCell *>solrefcell;//solution AMR cell
extern vector<OctCell *>solcoacell;//solution coase cell
extern set<OctCell*>ColCoa;

extern vector<Face> faces_comp;
extern Node *HeadListAllGrid;

static int N_total_adpcell=0;
static double Sigma_curl=0.0;
static double Sigma_div=0.0;//旋度和散度的标准差
static double Sigma_pio=0.0; //压力梯度
//计算自适应参数的均方根
//

extern Vec2D ExtPnt0;
extern Vec2D ExtPnt1;
extern Vec2D ExtPnt2;
extern Vec2D ExtPnt3;
extern double sumVorM;
//计算自适应参数 consider entropy

void getCellVorM(OctCell *unp)
{
    unp->vorM=(unp->Ux[2]-unp->Uy[1])*unp->dx*unp->dy;
}

void solveGradient();
void getAllVorM()
{
    solveGradient();

    Node *current;
    current = HeadListAllGrid;  
    OctCell *lsbl;

    Vec2D tpvec;
    double jl0,jl1,jl2,jl3;
    double maxhxy;
    maxhxy=max2(hx,hy);
    sumVorM=0.0;

    while(current != NULL)
    {		
        lsbl = current->cell;
        if(lsbl->flag == 0 && current->flg==0
            //   && lsbl->xc1<xOutlet-1.0*hx
            //    && lsbl->xc1>xInlet+1.0*hx
            //   && lsbl->yc1>ExtPnt1[1]+1.0*hy
            //   && lsbl->yc1<ExtPnt3[1]-1.0*hy
          ) {
#if IS_TUBE==1
            tpvec[0]=lsbl->xc1;
            tpvec[1]=lsbl->yc1;
            jl0=fabs(tpvec-ExtPnt0);
            jl1=fabs(tpvec-ExtPnt1);
            jl2=fabs(tpvec-ExtPnt2);
            jl3=fabs(tpvec-ExtPnt3);

            if(jl0>3.0*maxhxy && jl1>3.0*maxhxy
                && jl2 > 3.0*maxhxy && jl3>3.0*maxhxy){
#endif
                getCellVorM(lsbl);

                sumVorM+=lsbl->vorM;

#if IS_TUBE==1
            }
#endif
        }        
        current = current->next;
    }
}

void getAdaptiveParameters(OctCell *unp)
{
    double li;
    li=pow(sqrt(unp->dx*unp->dy),1.5);
    unp->tcur=fabs(unp->Ux[2]-unp->Uy[1])*li;
    unp->vorM=(unp->Ux[2]-unp->Uy[1])*unp->dx*unp->dy;
    unp->tdiv=fabs(unp->Ux[1]+unp->Uy[2])*li;
//0:den,1:press, 2:entropy,3:
#if ADAPTIVE_INDICATOR==0
    unp->pio=sqrt((unp->Ux[0]*unp->Ux[0])+(unp->Uy[0]*unp->Uy[0]))*li;
#elif ADAPTIVE_INDICATOR==1 
    unp->pio=sqrt((unp->Ux[3]*unp->Ux[3])+(unp->Uy[3]*unp->Uy[3]))*li;
#elif ADAPTIVE_INDICATOR==2
    double ason2, a2;//ason2:pressure
    ason2=GAM11*(unp->dof[3][0]-0.5*(unp->dof[1][0]*unp->dof[1][0]
   /unp->dof[0][0]+unp->dof[2][0]*unp->dof[2][0]/unp->dof[0][0]));
    a2=GAMMA*ason2/unp->dof[0][0];
    unp->pio=(fabs(unp->Ux[3]-a2*unp->Ux[0])
        +fabs(unp->Uy[3]-a2*unp->Uy[0]))*li;
#elif ADAPTIVE_INDICATOR==3
#else
#endif
}
extern double sumVorM;
void getAllAdpParametersList()
{
    Node *current;
    current = HeadListAllGrid;  
    OctCell *lsbl;

    Vec2D tpvec;
    double jl0,jl1,jl2,jl3;
    double maxhxy;
    Sigma_curl=0.0;
    Sigma_div=0.0;
    Sigma_pio=0.0;
    N_total_adpcell=0;
    maxhxy=max2(hx,hy);

    while(current != NULL)
    {		
        lsbl = current->cell;
        if(lsbl->flag == 0 && current->flg==0
            //   && lsbl->xc1<xOutlet-1.0*hx
            //    && lsbl->xc1>xInlet+1.0*hx
            //   && lsbl->yc1>ExtPnt1[1]+1.0*hy
            //   && lsbl->yc1<ExtPnt3[1]-1.0*hy
          ) {
#if IS_TUBE==1
            tpvec[0]=lsbl->xc1;
            tpvec[1]=lsbl->yc1;
            jl0=fabs(tpvec-ExtPnt0);
            jl1=fabs(tpvec-ExtPnt1);
            jl2=fabs(tpvec-ExtPnt2);
            jl3=fabs(tpvec-ExtPnt3);

            if(jl0>3.0*maxhxy && jl1>3.0*maxhxy
                && jl2 > 3.0*maxhxy && jl3>3.0*maxhxy){
#endif
                getAdaptiveParameters(lsbl);

                Sigma_curl += lsbl->tcur * lsbl->tcur;
                Sigma_div += lsbl->tdiv *lsbl->tdiv;			
                Sigma_pio+=lsbl->pio* lsbl->pio;

                ++N_total_adpcell;
#if IS_TUBE==1
            }
#endif
        }        
        current = current->next;
    }
    Sigma_curl = sqrt(Sigma_curl/N_total_adpcell);
    Sigma_div = sqrt(Sigma_div/N_total_adpcell);
    Sigma_pio = sqrt(Sigma_pio/N_total_adpcell);
}

//判断哪些网格需要加密
//int Namr=0;
void charg_sol_grid(OctCell *pg)//pg为网格单元的指针
{
	bool refflag;
    Vec2D tpvec;
    double jl0,jl1,jl2,jl3;
    double maxhxy;
    maxhxy=max2(hx,hy);
	switch(AMRP)
	{	
	case 0: 
        refflag=(pg->tcur > refine_coe*Sigma_curl
            || pg->tdiv > refine_coe1*Sigma_div);
        break;
	case 1: refflag=(pg->tdiv > refine_coe*Sigma_div); break;	
	case 2: refflag=(pg->tcur > refine_coe*Sigma_curl); break;
	case 3: refflag=(pg->pio > refine_coe*Sigma_pio); break;
	case 4: 
        refflag=(pg->tdiv > refine_coe*Sigma_div
            || pg->pio > refine_coe1*Sigma_pio);
        break;
    case 5:
        refflag=(pg->tdiv > refine_coe*Sigma_div || pg->trbamr);
        break;
	}

	if(refflag==true)
	{
#if IS_TUBE==1
        tpvec[0]=pg->xc1;
        tpvec[1]=pg->yc1;
        jl0=fabs(tpvec-ExtPnt0);
        jl1=fabs(tpvec-ExtPnt1);
        jl2=fabs(tpvec-ExtPnt2);
        jl3=fabs(tpvec-ExtPnt3);

        if(jl0>3.0*maxhxy && jl1>3.0*maxhxy
            && jl2 > 3.0*maxhxy && jl3>3.0*maxhxy)
#endif
            solrefcell.push_back(pg);

		//solrefcell[Namr]=pg;//用于存放需要加密的网格单元地址的数组
	    //Namr++;//记录需要加密的网格单元的数目
	}
}

void charg_sol_gridAll()//判断哪些网格单元需要加密
{
	Node *current;
	OctCell *lsbl;
	current = HeadListAllGrid;
//	Namr=0; //需要加密网格单元数目初始化
    solrefcell.clear();
	while(current != NULL)
	{
		lsbl = current->cell;
		if(current->flg==0 && lsbl->flag==0
#if NAR_MAXAMR==1
            && lsbl->level < lsbl->level0+Nar
#elif NAR_MAXAMR==0
            && lsbl->level < MaxAMR
#endif
            ) //最多只作MaxAMR次解自适应
            //20121111---------------------
  //          && lsbl->level < lsbl->level0+Nar) //最多只作Nar次解自适应
        {
            charg_sol_grid(lsbl);
        }
		current = current->next;
	}
}

//这个函数可以优化一下,对于pp使用set这样可以杜绝找到相同的地址值
//找到需要加密的网格单元的周边网格(一层)
void find_sol_refp(OctCell *pp[],OctCell *parent)
{ 		
	pp[0]=EastNeighbor(parent);	
	pp[1]=WestNeighbor(parent);		
	pp[2]=NorthNeighbor(parent);
	pp[3]=SouthNeighbor(parent);	
	pp[4]=WestNeighbor(pp[2]);
	pp[5]=EastNeighbor(pp[2]);		
	pp[6]=WestNeighbor(pp[3]); 	
	pp[7]=EastNeighbor(pp[3]);
}

void find_sol_4nbr(OctCell *pp[],OctCell *parent)
{ 		
	pp[0]=EastNeighbor(parent);	
	pp[1]=WestNeighbor(parent);		
	pp[2]=NorthNeighbor(parent);
	pp[3]=SouthNeighbor(parent);	
	//pp[4]=WestNeighbor(pp[2]);
	//pp[5]=EastNeighbor(pp[2]);		
	//pp[6]=WestNeighbor(pp[3]); 	
	//pp[7]=EastNeighbor(pp[3]);
}

int chargeleveldifference(OctCell *parent)
{
	int i;
	OctCell *pp[8]={NULL}; 
	int plevel,nlevel;//目标单元的层次-1，邻居单元的层次
	plevel=parent->level+1;
	find_sol_refp(pp,parent);
	for(i=0;i<8;i++)
	{
        if(pp[i]!=NULL)//here don't need if TODO
		{
			nlevel=pp[i]->level;
			if(abs(nlevel-plevel)>1)
			{break;}
		}
	}
	if(i==8) return(0);
	else return(1);
}

void charggrid(OctCell *pg);
void charggrid0(OctCell *pg);

void childrenvalue(OctCell *p1, OctCell *parent,int m,double hcx, double hcy);
void childrenvalue(OctCell *p1,OctCell *parent,double hx,double hy,int m);
void son_cell_value(OctCell *pcell, int nson, OctCell *pchild);
//由父亲pcell得到它的孩子的moment使用L2投影
void parent_cell_value(OctCell *pchildren[4], OctCell *pcell);
//由四个孩子得到父亲的moment
void creat_sol_children(OctCell *parent)
{
	int j,m,preflag,dlevel;
	OctCell *children[4]={NULL};
	//double plevel,hx,hy;
	double hx1,hy1;
	if(parent!=NULL)//here do not need if TODO
	{
		preflag=parent->reflag;	
		dlevel=chargeleveldifference(parent);
		if(preflag==1||dlevel==1) return;//判断该网格单元是否加密
		else
		{
		//	plevel=parent->level+1;		
		//	hx=Hx/pow(2.0,plevel);hy=Hy/pow(2.0,plevel);//子网格的步长	
            hx1=0.5*parent->dx;
            hy1=0.5*parent->dy;
			for(m=0;m<4;++m)
			{
				//children[m]=(struct cell*)malloc(LENC);
                children[m]=new OctCell(hx1,hy1);

				parent->children[m]=children[m];
				childrenvalue(children[m],parent,hx1,hy1,m);
				charggrid0(children[m]);// here I omit is right TODO
			}
			parent->reflag=1;
//            parent->p_node=NULL;
		}
	}

	for(m=0;m<4;m++)
	{
		children[m]->level0=parent->level0;//
        son_cell_value(parent,m,children[m]);
	}				
}
//加密网格
void creat_amr_grid()
{   
    int kr,Nk;
    int i,j,k,ir;
    OctCell *parent;

    // Nk=Namr-1;//需要加密的网格
    Nk=solrefcell.size()-1;
    //printf("refine solution grid number: Nk=%d\n",Nk);		   

    OctCell *pp[9]={NULL},*pc[9]={NULL},*p=NULL;
    //物面相交的网格邻近的需要加密的网格单元地址 
#if NAR_MAXAMR==1
    for(ir=0;ir<=Nr+Nar;ir++)
#elif NAR_MAXAMR==0
    for(ir=0;ir<MaxAMR;ir++)
#endif
    { 
        for(kr=0;kr<=Nk;kr++)//Nk与解自适应网格数	 		   
        {			   
            parent = solrefcell[kr];			  	 
          //  if(parent==NULL){cout<<"parent NULL\n";}
            if(parent->level==ir)
            {		   		   
                //找到需加密的网格邻近的需要加密的网格单元		  
                find_sol_refp(pp,parent);
               // find_sol_4nbr(pp,parent);
                pp[8]=parent;	 
                k=0;	   
                for(i=0;i<9;i++)		
                {
                    if(pp[i]!=NULL)//here do not need if TODO 
                    {pc[k]=pp[i];k++;}	  
                }
                //根据网格单元的层次对这些网格进行排序  	 
                for(j=1;j<=k-1;j++)	 
                {		  
                    for(i=1;i<=k-j;i++)				 
                        if(pc[i-1]->level>pc[i]->level)			   
                        {p=pc[i-1];pc[i-1]=pc[i];pc[i]=p;}
                }
                for(i=0;i<k;i++)
                {
                    if(pc[i]!=NULL&&pc[i]->reflag == 0)
                    {		 			  				 
                        //creat_sol_children(pc[i],pc[i]->cflag); 
                        creat_sol_children(pc[i]); 
                        //					   pc[i]->reflag=1;			  			  			   			
                    }	  		 		   	 
                }	
            }
        }
    }
    cout<<"refine solution grid number: Nk= "<<Nk<<endl;		   
}

/*
void creat_amr_grid(int iccc)
{   
    int kr,Nk;
    int i,j,ir;
    OctCell *parent;

    // Nk=Namr-1;//需要加密的网格
    Nk=solrefcell.size()-1;
    //printf("refine solution grid number: Nk=%d\n",Nk);		   

    cout<<"refine solution grid number: Nk= "<<Nk<<endl;		   
    const int  k=9;
    for(kr=0;kr<=Nk;kr++)//Nk与解自适应网格数	 		   
    {			   
        OctCell *pp[9]={NULL},*p=NULL;
        //物面相交的网格邻近的需要加密的网格单元地址 
        parent = solrefcell[kr];			  	 
        for(ir=0;ir<=Nr+Nar;ir++)
        { 
            if(parent->level==ir)
            {		   		   
                //找到需加密的网格邻近的需要加密的网格单元		  
                find_sol_refp(pp,parent);
                pp[8]=parent;	 
                //根据网格单元的层次对这些网格进行排序  	 
                //
                for(j=1;j<=k-1;j++)	 
                {		  
                    for(i=1;i<=k-j;i++)				 
                        if(pp[i-1]->level>pp[i]->level)			   
                        {p=pp[i-1];pp[i-1]=pp[i];pp[i]=p;}
                }
                for(i=0;i<k;i++)
                {		 
                    if(pp[i]!=NULL&&pp[i]->reflag == 0)	
                //    if(pp[i]->reflag == 0)	 
                    {		 			  				 
                        //creat_sol_children(pc[i],pc[i]->cflag); 
                        creat_sol_children(pp[i]); 
                        //					   pc[i]->reflag=1;			  			  			   			
                    }	  		 		   	 
                }	
                break;
            }
        }
    }

    cout<<"creat_amr_grid(int iccc)"<<endl;
}
*/

//判断哪些网格需要粗化
//int Ncoar=0;
void chargsolgridforcoarse(OctCell *pg)//pg为网格单元的指针
{
    //	int refflag;
    bool refflag;
    switch(COAR)
    {
    case 0: refflag=(pg->tcur <coarse_coe*Sigma_curl
                && pg->tdiv <coarse_coe1*Sigma_div);
            break;
    case 1: refflag=(pg->tdiv <coarse_coe*Sigma_div); break;
    case 2: refflag=(pg->tcur <coarse_coe*Sigma_curl); break;
    case 3: refflag=(pg->pio <coarse_coe*Sigma_pio); break;
    case 4: refflag=(pg->tdiv <coarse_coe*Sigma_div
                && pg->pio <coarse_coe1*Sigma_pio);
            break;
    case 5:
            refflag=((!pg->trbamr)
                && pg->tdiv <coarse_coe*Sigma_div);
            if(pg->exflag!=0) refflag=true;
            break;
    }
#if IS_TUBE==1
//    if(pg->exflag!=0) refflag=true;//20121112
#endif
    if(refflag==true && pg->level>pg->level0)
    {
       // solcoacell[Ncoar]=pg;//用于存放需要粗化的网格单元地址的数组
        solcoacell.push_back(pg);
        pg->coflag=true;
       // Ncoar++;//记录需要粗化的网格单元的数目
    }
}

void chargegridforcoarse()
{
	Node *current;
	OctCell *lsbl;
	current = HeadListAllGrid;
    solcoacell.clear();
	//Ncoar=0; //需要粗化网格单元数目初始化
	while(current != NULL)
	{
		lsbl = current->cell;
#if IS_TUBE==0
		if(lsbl->flag== 0
#elif IS_TUBE==1
        if((lsbl->flag== 0 || lsbl->exflag!=0)
#endif
            && current->flg<=2)
			chargsolgridforcoarse(lsbl);
		current = current->next;
	}
}

//删除子网格单元，给父网格单元赋值
void creatsolcoarseparent(OctCell *parent,OctCell *pp[4])
{
    for(int i=0;i<4;++i) {
        parent->dof[i][0]=0.25*(pp[0]->dof[i][0]+pp[1]->dof[i][0]
            +pp[2]->dof[i][0]+pp[3]->dof[i][0]);
      
        parent->dof[i][1]=0.125*(pp[0]->dof[i][1]+pp[1]->dof[i][1]
            +pp[2]->dof[i][1]+pp[3]->dof[i][1])
            +0.375*(-pp[0]->dof[i][0]+pp[1]->dof[i][0]
                -pp[3]->dof[i][0]+pp[2]->dof[i][0]);

        parent->dof[i][2]=0.125*(pp[0]->dof[i][2]+pp[1]->dof[i][2]
            +pp[2]->dof[i][2]+pp[3]->dof[i][2])
            +0.375*(-pp[0]->dof[i][0]-pp[1]->dof[i][0]
                +pp[3]->dof[i][0]+pp[2]->dof[i][0]);

        parent->dof[i][3]=0.0625*(pp[0]->dof[i][3]+pp[1]->dof[i][3]
            +pp[2]->dof[i][3]+pp[3]->dof[i][3])
            +0.1875*(-pp[0]->dof[i][2]+pp[1]->dof[i][2]
                -pp[3]->dof[i][2]+pp[2]->dof[i][2])
            +0.1875*(-pp[0]->dof[i][1]-pp[1]->dof[i][1]
                +pp[3]->dof[i][1]+pp[2]->dof[i][1])
            +0.5625*(pp[0]->dof[i][0]-pp[1]->dof[i][0]
                -pp[3]->dof[i][0]+pp[2]->dof[i][0]);

        parent->dof[i][4]=0.0625*(pp[0]->dof[i][4]+pp[1]->dof[i][4]
            +pp[2]->dof[i][4]+pp[3]->dof[i][4])
            +0.46875*(-pp[0]->dof[i][1]+pp[1]->dof[i][1]
                -pp[3]->dof[i][1]+pp[2]->dof[i][1]);

        parent->dof[i][5]=0.0625*(pp[0]->dof[i][5]+pp[1]->dof[i][5]
            +pp[2]->dof[i][5]+pp[3]->dof[i][5])
            +0.46875*(-pp[0]->dof[i][2]-pp[1]->dof[i][2]
                +pp[3]->dof[i][2]+pp[2]->dof[i][2]);
    }
    parent->level0=pp[0]->level0;
    for(int i=0; i<4;++i){
        parent->children[i]=NULL;
      if(pp[i]!=NULL)delete pp[i];//TODO
//        if(pp[i]!=NULL)ColCoa.insert(pp[i]);//TODO
    }
}
/*void release_ColCoa()
{
    for(set<OctCell*>::iterator iz=ColCoa.begin();
        iz!= ColCoa.end();++iz){
        delete *iz;
    }
    ColCoa.clear();
}
*/

int chargecoarselevel(OctCell *pc)
{
	int i;
	OctCell *pp[8]={NULL};
	find_sol_refp(pp,pc);
	for(i=0;i<4;i++)
	{
		if(pp[i]!=NULL&&pp[i]->reflag==1) 
		//if(pp[i]->reflag==1) 
		{
			break;
		}
	}
	if(i==4) return(0);
	else return(1);
}

//粗化网格
void creatcoarsegrid()
{   
    int kr,Nk;
    int i,ppflag;
    OctCell *parent,*p;
   // Nk=Ncoar-1;//需要加密的网格
   Nk=solcoacell.size()-1;
    //printf("coarse solution grid number: Nk=%d\n",Nk);		   

    for(kr=0;kr<=Nk;kr++)//Nk与解自适应网格数	 		  
    {		 		   		  		  
        OctCell *pp[4];
        //需要粗化的网格单元的兄弟网格单元的地址网格单元地址        		  			   
        p= solcoacell[kr];
        parent=p->parent;
        //	   printf("Kr=%d\n",kr);
        ppflag=0;
        if(parent!=NULL&&parent->reflag==1)
       // if(parent->reflag==1)
        {	   
          //  NOchild=p->NOchildren;
            for(i=0;i<4;i++)
            {
                pp[i]=parent->children[i];
                ppflag+=chargecoarselevel(pp[i]);
            } 
            if(ppflag==0 && pp[0]->coflag==true 
                && pp[1]->coflag==true && pp[2]->coflag==true 
                && pp[3]->coflag==true )//四个子节点都要粗华则粗化	   	   
            {				   
                creatsolcoarseparent(parent,pp);
                parent->coflag=false;
                parent->reflag=0;
            }
        }
    }
    cout<<"coarse solution grid number: Nk= "<<Nk<<endl;
}



/****************************
  判断是否是hole grid
 ****************************/
int chargeholegrid(OctCell *pc)
{
	OctCell *pp[8];
	find_sol_refp(pp,pc);
     if((pp[0]!=NULL&&pp[0]->reflag==1&&pp[1]!=NULL&&pp[1]->reflag==1)
     ||(pp[2]!=NULL&&pp[2]->reflag==1&&pp[3]!=NULL&&pp[3]->reflag==1))
    //if((pp[0]->reflag==1&&pp[1]->reflag==1) 
      //  ||(pp[2]->reflag==1&&pp[3]->reflag==1))
        return(1);
	else
		return(0);
}
/****************************
  判断是否是hill grid
 ****************************/
int chargehillgrid(OctCell *pc)
{
    OctCell *pp[8],*parent;
    parent=pc->parent;
    if(parent!=NULL)
    {
        if(parent->reflag==1)
        {
            find_sol_refp(pp,parent);
            if((pp[0]!=NULL&&pp[0]->reflag==0
                    &&pp[1]!=NULL&&pp[1]->reflag==0)
                &&(pp[2]!=NULL&&pp[2]->reflag==0
                    &&pp[3]!=NULL&&pp[3]->reflag==0))
     // if((pp[0]->reflag==0&&pp[1]->reflag==0)
      //          &&(pp[2]->reflag==0&&pp[3]->reflag==0))
                return(1);	
            else		
                return(0);
        }
    }
    return(0);
}

void formListForAllGrid();
void releaseList();
void setFlagForAll();
void setFlagForAll0();
//光顺网格
void smoothcell()
{
    //	long Nac;
    int hoflag,hiflag,i,nncell;
    Node *current;//新生成的四个节点
    OctCell *pc,*parent,*pp[4];//定义4个叶子节点的地址   
    int n_smth=3;
    do
    {
        //		releaselist(0);
        //		releasecelllist(Adpcelllist);
        //创建新的链表				
        //		Nac=FormAdpCellList(cartcell,cNx,cNy);

        releaseList();
        formListForAllGrid();
     //   setFlagForAll0();
        nncell=0;
        //current=Adpcelllist;
        current = HeadListAllGrid;

        while(current!=NULL)
        {             			
            pc=current->cell;
            if(pc->flag==0 && current->flg==0){
            //if(current->flg!=6){
                hiflag=0;
                hoflag=0;
                hoflag=chargeholegrid(pc);//判断是否是hole grid
//cout<<"smooth\n";
                //if(hoflag==1){printf("hoflag=%d\n",hoflag);getchar();}
                hiflag=chargehillgrid(pc);//判断是否是hill grid
                //if(hiflag==1){printf("hiflag=%d\n",hiflag);getchar();}
                /**/			
                if(hoflag==1) 
                {
                    creat_sol_children(pc);
                  //  pc->reflag=1;
                    nncell++;
                }

                if(hiflag==1)
                {	
                    parent=pc->parent;
                    for(i=0;i<4;i++)
                    {
                        pp[i]=parent->children[i];
                        pp[i]->coflag=true;
                    }
                    creatsolcoarseparent(parent,pp);
                    parent->coflag=false;
                    parent->reflag=0;
                    nncell++;
                }
            }
            current=current->next;
        }
        //	printf("nncell=%d\n",nncell);
        cout<<"nncell= "<<nncell<<endl;
        --n_smth;
    }while(nncell>0&& n_smth>0);
}

void solveEdge4Gradient(Face &fc, double fjbx[],double fjby[])
{
    OctCell *lsbl,*lsbl0;

    lsbl0=fc.parent;
    lsbl=fc.neighbor;
    double ush0[4],ush[4];
    double ujb0[4],ujb[4];
    for(int i=0;i<4;++i)
    {
        ush0[i]=lsbl0->dof[i][0];
        ush[i]=lsbl->dof[i][0];
    }
    shbl2jbbl(ush0,ujb0);
    shbl2jbbl(ush,ujb);

    for(int i=0; i<4;++i){
        fjbx[i]=0.5*(ujb0[i]+ujb[i])*fc.nml[0]*fc.area;
        fjby[i]=0.5*(ujb0[i]+ujb[i])*fc.nml[1]*fc.area;
    }
}

//void solveEdge4Gradient(Face &fc, double fjbx[],double fjby[])
void solveGradient()
{
    double fjbx[4],fjby[4];
    double vol0,vol1;

    for(vector<Face>::size_type isz=0; isz!=faces_comp.size(); ++isz){
        for(int i=0;i<4;++i){
            faces_comp[isz].parent->Ux[i]=0.0;
            faces_comp[isz].parent->Uy[i]=0.0;
            faces_comp[isz].neighbor->Ux[i]=0.0;
            faces_comp[isz].neighbor->Uy[i]=0.0;
        }
    }

    for(vector<Face>::size_type isz=0; isz!=faces_comp.size(); ++isz){
        solveEdge4Gradient(faces_comp[isz], fjbx,fjby);
        vol0=faces_comp[isz].parent->dx*faces_comp[isz].parent->dy;
        vol1=faces_comp[isz].neighbor->dx*faces_comp[isz].neighbor->dy;
        for(int i=0;i<4;++i){
            faces_comp[isz].parent->Ux[i]+=fjbx[i]/vol0;
            faces_comp[isz].parent->Uy[i]+=fjby[i]/vol0;;
            faces_comp[isz].neighbor->Ux[i]-=fjbx[i]/vol1;
            faces_comp[isz].neighbor->Uy[i]-=fjby[i]/vol1;
        }
    }
}

void son_cell_value(OctCell *pp, int nson, OctCell *pson)
{
    switch(nson) {
    case 0:
        for(int i=0;i<4;++i){
            pson->dof[i][0]=pp->dof[i][0]-0.5*pp->dof[i][1]
                -0.5*pp->dof[i][2]+0.25*pp->dof[i][3];

            pson->dof[i][1]=0.5*pp->dof[i][1]-0.25*pp->dof[i][3]
                -0.5*pp->dof[i][4];


            pson->dof[i][2]=0.5*pp->dof[i][2]-0.25*pp->dof[i][3]
                -0.5*pp->dof[i][5];

            pson->dof[i][3]=0.25*pp->dof[i][3];

            pson->dof[i][4]=0.25*pp->dof[i][4];
            pson->dof[i][5]=0.25*pp->dof[i][5];
        }
        break;
    case 1:
        for(int i=0;i<4;++i){
            pson->dof[i][0]=pp->dof[i][0]+0.5*pp->dof[i][1]
                -0.5*pp->dof[i][2]-0.25*pp->dof[i][3];

            pson->dof[i][1]=0.5*pp->dof[i][1]-0.25*pp->dof[i][3]
                +0.5*pp->dof[i][4];


            pson->dof[i][2]=0.5*pp->dof[i][2]+0.25*pp->dof[i][3]
                -0.5*pp->dof[i][5];

            pson->dof[i][3]=0.25*pp->dof[i][3];

            pson->dof[i][4]=0.25*pp->dof[i][4];
            pson->dof[i][5]=0.25*pp->dof[i][5];
        }
        break;
    case 2:
        for(int i=0;i<4;++i){
            pson->dof[i][0]=pp->dof[i][0]+0.5*pp->dof[i][1]
                +0.5*pp->dof[i][2]+0.25*pp->dof[i][3];

            pson->dof[i][1]=0.5*pp->dof[i][1]+0.25*pp->dof[i][3]
                +0.5*pp->dof[i][4];


            pson->dof[i][2]=0.5*pp->dof[i][2]+0.25*pp->dof[i][3]
                +0.5*pp->dof[i][5];

            pson->dof[i][3]=0.25*pp->dof[i][3];

            pson->dof[i][4]=0.25*pp->dof[i][4];
            pson->dof[i][5]=0.25*pp->dof[i][5];
        }
        break;
    case 3:
        for(int i=0;i<4;++i){
            pson->dof[i][0]=pp->dof[i][0]-0.5*pp->dof[i][1]
                +0.5*pp->dof[i][2]-0.25*pp->dof[i][3];

            pson->dof[i][1]=0.5*pp->dof[i][1]+0.25*pp->dof[i][3]
                -0.5*pp->dof[i][4];


            pson->dof[i][2]=0.5*pp->dof[i][2]-0.25*pp->dof[i][3]
                +0.5*pp->dof[i][5];

            pson->dof[i][3]=0.25*pp->dof[i][3];

            pson->dof[i][4]=0.25*pp->dof[i][4];
            pson->dof[i][5]=0.25*pp->dof[i][5];
        }
        break;
    }
}


void outputcell_sol(int ij=0, const string& str123="amr");
void outputcell_new(int ij=0);
void setFlagForAll();
void formListFaceCellAll();   
void face_vec();
void point_vec();     
void solveGradient();

extern PXYZ ext_wall[NPT_EXT_WALL+1]; //points at external computational domain
void formExtWallInOutBndry();
void set_exFlag4All(int nmpnt, PXYZ wallpt[]);//网格中心在固壁内部的为-1流场内部为0 在翼型上的为2, 尖后缘为4
void find_multivalue_cell(const PXYZ &sharppoint);
//void find_multivalue_cell(PXYZ &sharppoint);
void releaseFaceCellList();
extern PXYZ MBsharp_point[][8];
void kxrcf_relate();
void  face_kxrcf_all();
void  set_kxrcf();
void release_ColCoa();
void SolutionAMR(int i)
{	
#ifndef KXRCF_AMR
    solveGradient();
	getAllAdpParametersList();
#else
    kxrcf_relate();
    face_kxrcf_all();
    set_kxrcf();
#endif
	charg_sol_gridAll();
	chargegridforcoarse();
	creatcoarsegrid();
	creat_amr_grid();
//    release_ColCoa();
    solrefcell.clear();//solution AMR cell
    solcoacell.clear();//solution coase cell

	cout<<"Adaptively Refine and Coarsen Mesh Complete>>>>>\n";
	smoothcell();		
	cout<<"Smooth Adaptively Mesh Complete>>>>>\n";
    releaseList();
    formListForAllGrid();

    setFlagForAll();

    releaseFaceCellList();
    formListFaceCellAll();   
#if IS_TUBE==1
    set_exFlag4All(NPT_EXT_WALL, ext_wall);
    formExtWallInOutBndry();
#endif
    face_vec();
    point_vec();     
 // outputcell_new(i);//output mesh
  //outputcell_sol(i);//output solution
#if HAVE_MV!=0
    for(int wy=0; wy<N_BODY;++wy){
        for(int i_t=0;i_t<MVptNum[wy];++i_t){  
            find_multivalue_cell(MBsharp_point[wy][i_t]);
        }
    }
#endif
}
