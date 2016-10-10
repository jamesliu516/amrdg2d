
#include"non_uniform_grid.h"
//extern OctCell *bodygrid;
extern Node *HeadListAllGrid;
extern double dt;
extern double timesum;

void cal_te_res()
{
    Node *current;
    current = HeadListAllGrid;  
    OctCell *pcell0;
    double db_tmp=0.0;
    while(current != NULL)
    {		
        pcell0 = current->cell;
        if(pcell0->flag == 0 && current->flg<=2 ) {
            db_tmp+=pcell0->residual[0][0]*pcell0->residual[0][0]
                *pcell0->dx*pcell0->dy;
        }        
        current = current->next;
    }
    db_tmp=sqrt(db_tmp);
//    cout<<"rhoE ie. residual[3][0]->"<<scientific<<db_tmp<<"\n";
    cout<<"rho       "<<scientific<<db_tmp<<"\n";
}

