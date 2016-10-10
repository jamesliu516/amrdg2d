
#include"non_uniform_grid.h"
#include"findneighbor.h"

int newFlagpoint(double x, double y, int NumberAirfoilPoint,PXYZ wallpoint[] );
    //应用射线方法 //in the solid flag==-1, out the solid flag=0
    //at the surface of the solid flag=2
extern Node *HeadListAllGrid;
extern PXYZ MultiCircle[N_BODY+1][500]; //points at solid wall
extern int NWallPts[N_BODY+3];
//int newFlagpoint(double x, double y);////应用射线方法
//网格中心在airfoil中的为-1流场内部为0 在翼型上的为2, 尖后缘为4
//网格中心在流体中的为0,网格中心在固体内部且和流体网格相邻的为2其他网格为-1, TODO 多值点
void setFlagForAll()
{
    Node *current;
    OctCell *lsbl;
    OctCell *lss;
    OctCell *lse;
    OctCell *lsw;
    OctCell *lsn;
    double x1,y1;
    int value, flg[N_BODY+1];

    current = HeadListAllGrid;

    double xcn[4],ycn[4];
    int flgn[4], sflg;	
    int i;
    int n0,n2,nm;

    while(current != NULL)
    {
        lsbl = current->cell;

        x1=lsbl->xc1;
        y1=lsbl->yc1;

        for (int jj=0; jj<N_BODY;++jj){
            value = newFlagpoint(x1,y1,NWallPts[jj], MultiCircle[jj]);
            if(value==0) {
                flg[jj]=value;
            }
            else if(value==-1) {
                lss=SouthNeighbor(lsbl);
                lse=EastNeighbor(lsbl);
                lsw=WestNeighbor(lsbl);
                lsn=NorthNeighbor(lsbl);
                xcn[0]=lss->xc1;
                ycn[0]=lss->yc1;
                xcn[1]=lse->xc1;
                ycn[1]=lse->yc1;	
                xcn[2]=lsn->xc1;
                ycn[2]=lsn->yc1;
                xcn[3]=lsw->xc1;
                ycn[3]=lsw->yc1;	        
                sflg=0;
                for(i=0; i<4; ++i) {      
                    flgn[i]=newFlagpoint(xcn[i],ycn[i],NWallPts[jj], MultiCircle[jj]);

                    if(flgn[i]==0) {
                        flg[jj]=2;
                        ++sflg;
                        break;
                    }
                }  

                if(sflg==0) flg[jj]=-1;       
            }		       
        }
        n0=n2=nm=0;

        for (int jj=0; jj<N_BODY;++jj){
            if(flg[jj]==0) ++n0;
            if(flg[jj]==2) ++n2;
            if(flg[jj]==-1) ++nm;
        }

        if(n0==N_BODY) lsbl->flag=0;
        else if(n2>0) lsbl->flag=2;
        else lsbl->flag=-1;
        current=current->next;
    }
}
//
//网格中心在流体中的为0,网格中心在固体内部且和流体网格相邻的为2其他网格为-1, TODO 多值点
//in order to solve internal flow problem
////网格中心在固壁内部的为-1流场内部为0 在翼型上的为2, 尖后缘为4
void set_exFlag4All(int nmpnt, PXYZ wallpt[]) 
{
    Node *current;
    OctCell *lsbl;
    OctCell *lsbl1[4];
    double x1,y1;
    int value, flg=0;

    double xcn[4],ycn[4];
    int flgn[4], sflg;	

    current = HeadListAllGrid;
    while(current != NULL)
    {
        lsbl = current->cell;

        x1=lsbl->xc1;
        y1=lsbl->yc1;
        flg=0;

        value = newFlagpoint(x1,y1,nmpnt, wallpt);
        if(value==0 ) {
            flg=-1;
            lsbl->flag=-1; ///
        }
        lsbl->exflag=flg;
        current=current->next;
    }
    //cout<<"it is ok\n";

    current = HeadListAllGrid;
    while(current != NULL)
    {
        lsbl = current->cell;
        if(lsbl->exflag==-1){
            if(current->flg!=6)
            {
                //-----------------------
                lsbl1[0]=SouthNeighbor(lsbl);
                lsbl1[1]=EastNeighbor(lsbl);	
                lsbl1[2]=NorthNeighbor(lsbl);
                lsbl1[3]=WestNeighbor(lsbl);	        
                sflg=0;
                for(int i=0; i<4; ++i) {      
                    if(lsbl1[i]->exflag==0) {
                        flg=2;
                        ++sflg;
                        break;
                    }
                }  

                if(sflg==0) flg=-1;      
                lsbl->exflag=flg;
            }
        }
        current=current->next;
    }
}

void charggrid0(OctCell*);
//flag=0:fluid cell, flag=-1:solid cell, flag=2:interface cell(只要相交)
void setFlagForAll0() 
{
    Node *current;
    OctCell *lsbl;
    double x1,y1;
    int value, flg;

    current = HeadListAllGrid;

    double xcn[4],ycn[4];
    int flgn[4], sflg;	

    while(current != NULL)
    {
        lsbl = current->cell;

        charggrid0(lsbl);
      //  lsbl->flag=flg; //here wrong 2012.4.20
        current=current->next;
    }
    cout<<"setFlagAll0 finished\n";
}

