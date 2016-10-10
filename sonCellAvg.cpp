//get the cell 4 son's average of a cell
#include"non_uniform_grid.h"
void sonCellAvg(int nson, OctCell *pcell, double dof0[4])
{
    double dof[4][6]={0.0};
    for(int i=0;i<4;++i){
        for(int j=0;j<nDOF;++j){
            dof[i][j]=pcell->dof[i][j];
        }
    }

    switch(nson) {
    case 0:
        for(int i=0;i<4;++i) {
            dof0[i]= dof[i][0]-0.5*dof[i][1]-0.5*dof[i][2]+0.25*dof[i][3];
        }
        break;
    case 1:
        for(int i=0;i<4;++i){
            dof0[i]= dof[i][0]+0.5*dof[i][1]-0.5*dof[i][2]-0.25*dof[i][3];
        }
        break;
    case 2:
        for(int i=0;i<4;++i){
            dof0[i]= dof[i][0]+0.5*dof[i][1]+0.5*dof[i][2]+0.25*dof[i][3];
        }
        break;
    case 3:
        for(int i=0;i<4;++i){
            dof0[i]= dof[i][0]-0.5*dof[i][1]+0.5*dof[i][2]-0.25*dof[i][3];
        }
        break;
    }
}

