//由叶子节点获得他的父亲的单元平均,好用于限制器设置,或者将来的多重网格
#include"non_uniform_grid.h"
#include"findneighbor.h"

extern Node *HeadListAllGrid;
//pcell 是叶子节点
const int back_para=4; //最大设置平均到爷爷的父辈

void limit_cell_avg(OctCell *pcell);

void setAllLeafbAvgFalse()
{
    Node *current;
    current = HeadListAllGrid;  
    OctCell *pcell0;

    current = HeadListAllGrid;     
    while(current != NULL)
    {		
        pcell0 = current->cell;
        if(pcell0->flag%2 == 0 && current->flg<=2 ) {
            pcell0->bAvg=false;
        }
        current = current->next;
    }			
}

void parent0bAvg(OctCell *pcell)
{
    int lev=pcell->level;
    OctCell *parent;
    parent=pcell->parent;
    int isum=0;

    // while(lev>=0) { //>0 or >=0 ?? TODO
    while(parent!=NULL) {
        if(parent->bAvg!=false){
            if(parent->children[0]->bAvg== false 
                && parent->children[1]->bAvg== false
                && parent->children[2]->bAvg== false
                && parent->children[3]->bAvg== false){
                parent->bAvg=false;
            }
            else{
                break;
            }

        }   
        ++isum;
        if(isum>back_para) break;
        parent=parent->parent;
    }        
}


void parent_avg(OctCell *pcell)
{
    int lev=pcell->level;

    OctCell *parent;
    parent=pcell->parent;
    int isum=0;

    // while(lev>=0) { //>0 or >=0 ?? TODO
    while(parent!=NULL) {
        if(parent->bAvg!=true){
            if(parent->children[0]->bAvg== true 
                && parent->children[1]->bAvg== true
                && parent->children[2]->bAvg== true 
                && parent->children[3]->bAvg== true){
                limit_cell_avg(parent); //由儿子的单元平均获得父亲的单元平均
            }
            else{
                break;
            }

        }   
        ++isum;
        if(isum>back_para) break;
        parent=parent->parent;

    }        
}   

//obtain the cell avarage from his four sons' cell moment
void limit_cell_avg(OctCell *pcell)
{
    pcell->bAvg=true;

    for(int i=0; i<4;++i) {
        pcell->dof[i][0]=0.25*(pcell->children[0]->dof[i][0] 
            + pcell->children[1]->dof[i][0] 
            + pcell->children[2]->dof[i][0] +pcell->children[3]->dof[i][0] );
    }
}

void tree_0bAvg()
{
    Node *current;
    current = HeadListAllGrid;  
    OctCell *pcell0;
    short lev0;
    OctCell *pnbr[4]={NULL};
    short levnbr[4];

    current = HeadListAllGrid;     
    while(current != NULL)
    {		
        pcell0 = current->cell;
        if(pcell0->flag%2 == 0 && current->flg<=2 ) {
            lev0=pcell0->level;

            pnbr[0]=EastNeighbor(pcell0);
            pnbr[1]=SouthNeighbor(pcell0);
            pnbr[2]=WestNeighbor(pcell0);
            pnbr[3]=NorthNeighbor(pcell0);

            for(int i=0;i<4;++i){
                levnbr[i]=pnbr[i]->level;
              //  reflg0[i]=pnbr[i]->reflag;
                if(lev0>levnbr[i]){
                    parent0bAvg(pcell0);     
                    //for(int j=0;j<4;++j){
                    //    parent0bAvg(pnbr[j]);     
                    //}
                    break;
                }
            }
        }
        current = current->next;
    }			
}

/*
   OctCell *NorthNeighbor(OctCell *pp);
   OctCell *WestNeighbor(OctCell *pp);
   OctCell *EastNeighbor(OctCell *pp);
   OctCell *SouthNeighbor(OctCell *pp);
   */

void tree_avg()
{
    Node *current;
    current = HeadListAllGrid;  
    OctCell *pcell0;
    OctCell *pnbr[4]={NULL};
    short levnbr[4]={0};
   // short reflgnbr[4]={0};
    short lev0;
    //short reflg0;

    current = HeadListAllGrid;     
    while(current != NULL)
    {		
        pcell0 = current->cell;
        if(pcell0->flag%2 == 0 && current->flg<=2 ) {
            lev0=pcell0->level;
           // reflg0=pcell0->reflag;

            pnbr[0]=EastNeighbor(pcell0);
            pnbr[1]=SouthNeighbor(pcell0);
            pnbr[2]=WestNeighbor(pcell0);
            pnbr[3]=NorthNeighbor(pcell0);
            for(int i=0;i<4;++i){
                levnbr[i]=pnbr[i]->level;
              //  reflg0[i]=pnbr[i]->reflag;
                if(lev0>levnbr[i]){
                    parent_avg(pcell0);     
                    // for(int j=0;j<4;++j){
                   //   parent_avg(pnbr[j]);     
                     //}
                    break;
                }
            }
        }
        current = current->next;
    }			
}

