#include"non_uniform_grid.h"
#include"findneighbor.h"
#include<set>
#include<map>
#include<cstddef>
extern OctCell *bodygrid;

//extern PXYZ Circle[NPOINTCIRCLE+1];

//const PXYZ *const wallpoint = Circle;
//const int NumberAirfoilPoint=NPOINTCIRCLE;
//const int HalfNumberAirfoilPoint=NumberAirfoilPoint / 2;
//const int n_sharp_point=HalfNumberAirfoilPoint;

void find_sol_refp(OctCell *pp[],OctCell *parent);

//PXYZ sharp_point=wallpoint[n_sharp_point];

extern set<OctCell*>mvalue_set;
void insert_cell2(OctCell *pp, set<OctCell*> &grid)
{
    if(pp->reflag==0){
        if(pp->flag%2==0)grid.insert(pp);
    }
    else {
        for(int m=0; m<4;++m) 
            insert_cell2(pp->children[m], grid);
    }
}

//这里只是找到可能有多值点情况的
void find_multivalue_cell(const PXYZ &sharppoint)
{
    double hcx, hcy;
    OctCell *cellflow=NULL;
    OctCell *in_this_cell=NULL;
    OctCell *ppvec[9]={NULL};
   // PXYZ point;
   // point=wallpoint[n_sharp_point];
    mvalue_set.clear();
    for(int i=1;i<=Nx*Ny;++i){
        cellflow=bodygrid+i;
        hcx=cellflow->dx;
        hcy=cellflow->dy;

        if(sharppoint.x < cellflow->xc1 + 0.5 * hcx + ERRS
            && sharppoint.x > cellflow->xc1 - 0.5 * hcx - ERRS
            && sharppoint.y < cellflow->yc1 + 0.5 * hcy + ERRS
            && sharppoint.y > cellflow->yc1 - 0.5 * hcy - ERRS){
            in_this_cell=cellflow;
            break;
        }
    }

    find_sol_refp(ppvec,in_this_cell);
    ppvec[8]=in_this_cell;
    for(int i=0; i<9;++i)
        insert_cell2(ppvec[i],mvalue_set);

    for(set<OctCell*>::iterator it_set=mvalue_set.begin();
        it_set!=mvalue_set.end();++it_set) {
        (*it_set)->mvflg=true;
    }
}


