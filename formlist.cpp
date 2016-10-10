
#include"non_uniform_grid.h"

extern Node *HeadListAllGrid;
//extern Node *HeadList1, *HeadListAllGrid, *HeadListAdp;
extern OctCell *bodygrid;
extern map<OctCell *, cellEdgeBool> cellEbools;
extern BndryNode *HeadFaceCell;
extern vector<OctCell*> ExtWallBndry;
extern vector<OctCell*> InletBndry;
extern vector<OctCell*> OutletBndry;
void formListAll(OctCell *pcell, int flg);

void formListForAllGrid()
{
	int j,k;
	int flg;
    cellEbools.clear();
		for(j=1; j<=Ny; j++)
		{
			for(k=1; k<=Nx; k++)
			{
			    if(  (k>=adp_begin+1 && k<= Nx-adp_begin)&&(j>=adp_begin+1&& j<= Ny-adp_begin) ){
	                      flg=0;
	                     }
	                  else if( (k<=2||k>=Nx-1) || (j<=2||j>=Ny-1)  ) {
	                         flg=6;
	                    }
	                   else {
	                        flg=2;
	                    }
	                    
			    formListAll(&bodygrid[(j-1)*Nx + k], flg );
			}
		}
        cout<<"all cell list has been formed\n";
}

static Node *Prev3=NULL;

void formListAll(OctCell *pcell, int flg)  //TODO cut off the solid cell
{
    cellEdgeBool ceb1;

	if(pcell!=NULL&&pcell->reflag==0)
	{
		Node *current;
	//    current = (Node *)malloc(sizeof(Node));
        current=new Node;
	
		if(HeadListAllGrid == NULL)	
		{	HeadListAllGrid = current;
            current->prev=NULL;
		}
		else	
			Prev3->next = current;
			
	    current->next = NULL;
	    current->flg=flg;	    
	    current->cell = pcell;
	    //设置叶子节点的 bool bAvg
	    pcell->bAvg=true;
        pcell->mvflg=false;
  //      pcell->p_node=current;
	    //
			      
	    cellEbools.insert(make_pair(pcell,ceb1));

	    current->prev=Prev3;
	    Prev3=current;
	}
	else
	    for(int m=0;m<4;m++) formListAll(pcell->children[m], flg);
}

extern Vec2D ExtPnt0;
extern Vec2D ExtPnt1;
extern Vec2D ExtPnt2;
extern Vec2D ExtPnt3;
void formExtWallInOutBndry()
{
	Node *current, *pnext;
	OctCell *pcell;
    ExtWallBndry.clear();
    InletBndry.clear();
    OutletBndry.clear();
	current = HeadListAllGrid;
    while(current != NULL)
    {
        pnext = current->next;
        pcell=current->cell;

        if(pcell->exflag==2 && (pcell->xc1>xInlet
                && pcell->xc1<xOutlet) ){
            ExtWallBndry.push_back(pcell);
        }

        if(pcell->exflag==2 && pcell->xc1<xInlet
            &&( pcell->yc1<ExtPnt0[1] && pcell->yc1>ExtPnt1[1] ) ){
            InletBndry.push_back(pcell);
        }

        if(pcell->exflag==2 && pcell->xc1>xOutlet
          &&(pcell->yc1<ExtPnt3[1] && pcell->yc1>ExtPnt2[1])){
            OutletBndry.push_back(pcell);
        }
        current=pnext;
    }
    cout<<"exteral wall face cells have been formed\n";
}


static BndryNode *PrevFacecell=NULL;
static int sfcell=0;
void formListFaceCell(OctCell *pg, int listflg)
{
	if(pg->flag==listflg)
	{
	    BndryNode *current;
	    current = new BndryNode;
	
	    if(HeadFaceCell == NULL)	
            HeadFaceCell = current;
	    else	
            PrevFacecell->next = current;
		
	    current->next = NULL;	  
	      
	    current->cell = pg;
	      
	    current->prev=PrevFacecell;
	    
//	    pg->pnode=current;
	    PrevFacecell=current;
	    ++sfcell;
	}
	
}

void formListFaceCellAll()
{
	Node *current, *pnext;
	OctCell *pcell;
       sfcell=0;
	current = HeadListAllGrid;
	while(current != NULL)
	{
		pnext = current->next;
		pcell=current->cell;
		
		formListFaceCell(pcell, 2);
		current=pnext;
	}
    cout<<"internal face cells have been formed\n";
}

//确定解自适应前所有网格单元的层数level0，解自适应加密在level0的基础上最多作Nar次
//粗化操作不能低于level0
void initiallevel()
{
	Node *currentc;
	OctCell *pc;
	currentc=HeadListAllGrid;
	while(currentc!=NULL)
	{
		pc=currentc->cell;

		pc->level0=pc->level;

		currentc=currentc->next;
	}
}

void releaseList()
{
    Node *current, *pnext;

    current = HeadListAllGrid;
    while(current != NULL)
    {
        pnext = current->next;
        delete current;
        current=pnext;
    }
    cout<<"release all listAllGrid end.\n";

    HeadListAllGrid=NULL;
    Prev3=NULL;

}
void releaseFaceCellList()
{
    BndryNode *curr, *pne;      
    curr=HeadFaceCell;
    while(curr!=NULL) {
        pne=curr->next;
    //    curr->cell->pnode=NULL;
        delete curr;
        curr=pne;
    }
    HeadFaceCell=NULL;
    PrevFacecell=NULL;
    cout<<"total internal face cell: "<<sfcell<<endl;
    cout<<"releasing internal face cell done\n";
    sfcell=0;
}

/*
Node *Prev1=NULL;

void formList1(OctCell *pcell)
{
	int m;
	if(pcell->reflag==0)
	{
		Node *current;
	    current = (Node *)malloc(sizeof(Node));
	
		if(HeadList1 == NULL)	
			HeadList1 = current;
		else	
			Prev1->next = current;
		current->next = NULL;
	    current->cell = pcell;
	    Prev1=current;
	}
	else
		for(m=0;m<4;m++) formList1(pcell->children[m]);
}

//这里这些链表的形成可以使用所有的网格链表去除到不要的单元就能得到新的链表这样做更省时
//可以优化程序
void formListForGrid(void)
{
	int i,j,k;
	for (i=1; i<=1; i++)
	{
		for(j=3; j<=Ny-2; j++)
		{
			for(k=3; k<=Nx-2; k++)
			{
				formList1(&bodygrid[(i-1)*Nx*Ny + (j-1)*Nx + k]);
			}
		}
	}
}
*/

/*
Node *Prev2=NULL;
void formListAdp(OctCell *pcell)
{
	int m;
	if(pcell->reflag==0)
	{
		Node *current;
	    current = (Node *)malloc(sizeof(Node));
	
		if(HeadListAdp == NULL)	
			HeadListAdp = current;
		else	
			Prev2->next = current;

		current->next = NULL;
	    current->cell = pcell;

	    Prev2=current;
	}
	else
		for(m=0;m<4;m++) formListAdp(pcell->children[m]);

}




void formListForAdp(void)
{
	int i,j,k;
	for (i=1; i<=1; i++)

	{
		for(j=5; j<=Ny-5; j++)
		{
			for(k=5; k<=Nx-5; k++)

			{
				formListAdp(&bodygrid[(i-1)*Nx*Ny + (j-1)*Nx + k]);
			}
		}

	}
}

*/

/*
void releaseListAdp(void)
{
	Node *current, *pnext;

	current = HeadListAdp;
	while(current != NULL)
	{
		pnext = current->next;
		free(current);
		current=pnext;
	}
	printf("release all listAdpGrid end.\n");
	HeadListAdp=NULL;
}*/

