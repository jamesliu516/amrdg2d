

#include"non_uniform_grid.h"

extern OctCell *bodygrid;
//extern double initQ,initP,uFree,vFree,initOmg,initK;

void shblToJbbl(SHBL *, JBBL *);
void jbblToShbl(JBBL *, SHBL *);
//come from wu z n'book,2008.2.21修改边界条件的法线方向,吴的书中法线是指向求解域内部
void farFieldBoundary(int itt)
{
    int j,k;
    int i1, j1;
    double nx;
    double ny;
    SHBL sh;
    JBBL jb1, jb2;

    double speed, a, vbn;

    vbn=0.0;
    //-----------------------------------------------------------------------

    //----远场边界条件，位于和y坐标轴垂直的下部面----begin--	 
    nx=0.0;
    ny = 1.0;

    for(k=1; k <= Nx; k++)
    {
        //----下部那个面-------begin--
        sh.q = bodygrid[2*Nx +k].dof[0][0];
        sh.qu = bodygrid[2*Nx +k].dof[1][0];
        sh.qv = bodygrid[2*Nx +k].dof[2][0];
        sh.te = bodygrid[2*Nx +k].dof[3][0];

        bodygrid[Nx +k].dof[0][0] = sh.q;
        bodygrid[Nx +k].dof[1][0]= sh.qu;
        bodygrid[Nx +k].dof[2][0]= sh.qv;
        bodygrid[Nx +k].dof[3][0]= sh.te;
        for( i1=0; i1<4; ++i1) 
            for(j1=1; j1<nDOF; ++j1)
                bodygrid[Nx +k].dof[i1][j1]=0.0;


        //---虚拟网格外插------
        for(i1=0; i1<4; ++i1) {
            bodygrid[k].dof[i1][0] = 2.0*bodygrid[Nx+k].dof[i1][0]
                - bodygrid[2*Nx+k].dof[i1][0];
            for(j1=1; j1<nDOF;++j1)
                bodygrid[k].dof[i1][j1] = 0.0;
        }

        //----下部那个面-------end--
    }


    //////////////////////////////////////////////////////////////////////////////////---------------------2008.7.1

    //     for(j=1;j<=Ny;j++) {
    //       for(k=1;k<=Nx;k++) {	         				       		       
    //			bodygrid[(j-1)*Nx +k]

    //----远场边界条件，位于和y坐标轴垂直的上部面----begin--	 
    nx=0.0;
    ny = -1.0;


    for(k=1; k <= Nx; k++)
    {
        //----上部那个面-------begin--
        sh.q = bodygrid[(Ny-3)*Nx+k].dof[0][0];
        sh.qu = bodygrid[(Ny-3)*Nx+k].dof[1][0];
        sh.qv = bodygrid[(Ny-3)*Nx+k].dof[2][0];
        sh.te = bodygrid[(Ny-3)*Nx+k].dof[3][0];


        bodygrid[(Ny-2)*Nx+k].dof[0][0] = sh.q;
        bodygrid[(Ny-2)*Nx+k].dof[1][0]= sh.qu;
        bodygrid[(Ny-2)*Nx+k].dof[2][0] = sh.qv;
        bodygrid[(Ny-2)*Nx+k].dof[3][0] = sh.te;

        for(i1=0; i1<4; ++i1)
            for(j1=1; j1<nDOF; ++j1) 
                bodygrid[(Ny-2)*Nx+k].dof[i1][j1] = 0.0;

        //--虚拟网格外插------
        for (i1=0; i1<4; ++i1) {
            bodygrid[(Ny-1)*Nx+k].dof[i1][0] = 2.0*bodygrid[(Ny-2)*Nx+k].dof[i1][0]
                - bodygrid[(Ny-3)*Nx+k].dof[i1][0];
            for(j1=1; j1<nDOF; ++j1)
                bodygrid[(Ny-1)*Nx+k].dof[i1][j1] = 0.0;
        }
        //----------------------------------------------
        //----上部那个面-------end--
    }


    //-------------------右侧面---------0000000000000------begin---
    nx = -1.0;
    ny=0.0;

    for(j=1; j <= Ny; j++)
    {
        /*-----------begin--*/
        sh.q = bodygrid[(j-1)*Nx+Nx-2].dof[0][0];
        sh.qu = bodygrid[(j-1)*Nx+Nx-2].dof[1][0];
        sh.qv = bodygrid[(j-1)*Nx+Nx-2].dof[2][0];
        sh.te = bodygrid[(j-1)*Nx+Nx-2].dof[3][0];

        bodygrid[(j-1)*Nx+Nx-1].dof[0][0]= sh.q;
        bodygrid[(j-1)*Nx+Nx-1].dof[1][0] = sh.qu;
        bodygrid[(j-1)*Nx+Nx-1].dof[2][0] = sh.qv;
        bodygrid[(j-1)*Nx+Nx-1].dof[3][0] = sh.te;

        for(i1=0; i1<4; ++i1)
            for(j1=1; j1<nDOF; ++j1)
                bodygrid[(j-1)*Nx+Nx-1].dof[i1][j1]=0.0;

        //---虚拟网格外插------
        for(i1=0; i1<4; ++i1) {
            bodygrid[(j-1)*Nx+Nx].dof[i1][0]
                = 2.0*bodygrid[(j-1)*Nx+Nx-1].dof[i1][0] - bodygrid[(j-1)*Nx+Nx-2].dof[i1][0];
            for(j1=1; j1<nDOF; ++j1)
                bodygrid[(j-1)*Nx+Nx].dof[i1][j1]=0.0;	
        }	       

    }

    //----右侧面------------------end--------

    //-------------------左侧面----begin------
    nx = 1.0;
    ny=0.0;

    for(j=1; j <= Ny; j++)
    {
        //-----------begin--------
        sh.q = bodygrid[(j-1)*Nx+3].dof[0][0];
        sh.qu = bodygrid[(j-1)*Nx+3].dof[1][0];
        sh.qv = bodygrid[(j-1)*Nx+3].dof[2][0];
        sh.te = bodygrid[(j-1)*Nx+3].dof[3][0];


        bodygrid[(j-1)*Nx+2].dof[0][0] = sh.q;
        bodygrid[(j-1)*Nx+2].dof[1][0]= sh.qu;
        bodygrid[(j-1)*Nx+2].dof[2][0] = sh.qv;
        bodygrid[(j-1)*Nx+2].dof[3][0]= sh.te;
        for(i1=0; i1<4; ++i1) 
            for(j1=1; j1<nDOF; ++j1)
                bodygrid[(j-1)*Nx+2].dof[i1][j1]=0.0;

        //---虚拟网格外插------
        for(i1=0; i1<4; ++i1) {
            bodygrid[(j-1)*Nx+1].dof[i1][0]
                = 2.0*bodygrid[(j-1)*Nx+2].dof[i1][0] - bodygrid[(j-1)*Nx+3].dof[i1][0];
            for(j1=1; j1<nDOF; ++j1)
                bodygrid[(j-1)*Nx+1].dof[i1][j1]=0.0;
        }
    }

    //----左侧面-------------------------end----
}

void getGridIndex(OctCell *parent, const PXYZ *p_sm, ptrOctCell &pcell);
void getGridIndex2(OctCell *parent, const PXYZ *p_sm, ptrOctCell &pcell);
void findIntpltCell(OctCell * incell, const PXYZ &gstp, const PXYZ &smtricp, vector< OctCell * > &grid) ;
void chazhi4smtrictPshbl(vector< OctCell * > &grid, const PXYZ &smtrict, JBBL &jbu );
//in2out can only set 0 (inlet) or 1 (outlet).
void getInletOutletBndry_1(OctCell *ghostcell, int in2out)
{
    double nx,ny;
    double speed, a, vbn=0.0;
    PXYZ symtrcpt; 
    PXYZ ghostpt, plsA, plsB;
    int i1,j1;

    //Rieman problem界面量
    //    double dil,dir,ui,prei;
    JBBL jb1,jb2; 

    SHBL sh;

    //	OctCell *leafp;
    OctCell *pincell=NULL;
    vector< OctCell * > grid;


    ghostpt.x = ghostcell->xc1;
    ghostpt.y = ghostcell->yc1;	
    if(in2out==0){
        symtrcpt.x=2.0*xInlet-ghostpt.x;
        symtrcpt.y=ghostpt.y;
        nx=1.0;
        ny=0.0;
    }
    else {
        symtrcpt.x=2.0*xOutlet-ghostpt.x;
        symtrcpt.y=ghostpt.y;
        nx=-1.0;
        ny=0.0;
    }
    getGridIndex2(ghostcell, &symtrcpt, pincell);
    findIntpltCell(pincell, ghostpt, symtrcpt, grid);
    chazhi4smtrictPshbl(grid, symtrcpt, jb1); 

    a = sqrt(GAMMA * jb1.p / jb1.q);

    speed = jb1.u * nx + jb1.v * ny;

    if(speed > 0.0 && fabs(speed) >= a)  //inflow
    {
        jb2.q=initQ;
        jb2.p=initP;
        jb2.u=uFree;
        jb2.v=vFree;

        jbblToShbl(&jb2, &sh);

        ghostcell->dof[0][0] = sh.q;
        ghostcell->dof[1][0]= sh.qu;
        ghostcell->dof[2][0]= sh.qv;
        ghostcell->dof[3][0]= sh.te;
        for( i1=0; i1<4; ++i1) 
            for(j1=1; j1<nDOF; ++j1)
                ghostcell->dof[i1][j1]=0.0;

    }
    else if(speed > 0.0 && fabs(speed) < a) //inflow
    {

        jb2.p = 0.5 * (initP + jb1.p) + 0.5 * jb1.q * a * (nx * (uFree - jb1.u) + ny * (vFree - jb1.v));
        vbn = 0.5 / (jb1.q * a) *  (initP - jb1.p) + 0.5 * (nx * (uFree + jb1.u) + ny * (vFree + jb1.v));

        jb2.q = initQ * pow(jb2.p/initP, 1./GAMMA);

        jb2.u = uFree + (vbn - (uFree * nx + vFree * ny )) * nx;
        jb2.v = vFree + (vbn - (uFree * nx + vFree * ny )) * ny;

        jbblToShbl(&jb2, &sh);

        ghostcell->dof[0][0] = sh.q;
        ghostcell->dof[1][0]= sh.qu;
        ghostcell->dof[2][0]= sh.qv;
        ghostcell->dof[3][0]= sh.te;
        for(i1=0; i1<4; ++i1) 
            for(j1=1; j1<nDOF; ++j1)
                ghostcell->dof[i1][j1]=0.0;

    }
    else if(speed <= 0.0 && fabs(speed) >= a) //outflow
    {
        jb2.q=jb1.q;
        jb2.p=jb1.p;
        jb2.u=jb1.u;
        jb2.v=jb1.v;

        jbblToShbl(&jb2, &sh);

        ghostcell->dof[0][0] = sh.q;
        ghostcell->dof[1][0]= sh.qu;
        ghostcell->dof[2][0]= sh.qv;
        ghostcell->dof[3][0]= sh.te;
        for(i1=0; i1<4; ++i1) 
            for(j1=1; j1<nDOF; ++j1)
                ghostcell->dof[i1][j1]=0.0;

    }
    else if(speed <= 0.0 && fabs(speed) < a) //outflow
    {

        jb2.p = 0.5 * (initP + jb1.p) + 0.5 * jb1.q * a * (nx * (uFree - jb1.u) + ny * (vFree - jb1.v));
        vbn = 0.5 / (jb1.q * a) *  (initP - jb1.p) + 0.5 * (nx * (uFree + jb1.u) + ny * (vFree + jb1.v));
        jb2.q = jb1.q * pow(jb2.p/jb1.p, 1./GAMMA);

        jb2.u = jb1.u + (vbn - (jb1.u * nx + jb1.v * ny)) * nx;

        jb2.v = jb1.v + (vbn - (jb1.u * nx + jb1.v * ny)) * ny;


        jbblToShbl(&jb2, &sh);

        ghostcell->dof[0][0] = sh.q;
        ghostcell->dof[1][0]= sh.qu;
        ghostcell->dof[2][0]= sh.qv;
        ghostcell->dof[3][0]= sh.te;
        for(i1=0; i1<4; ++i1) 
            for( j1=1; j1<nDOF; ++j1)
                ghostcell->dof[i1][j1]=0.0;

    }
}
extern JBBL State_l;
void getInletOutletBndry(OctCell *ghostcell, int in2out)
{
    double nx,ny;
    double speed, a, vbn=0.0;
    PXYZ symtrcpt; 
    PXYZ ghostpt, plsA, plsB;
    int i1,j1;

    //Rieman problem界面量
    //    double dil,dir,ui,prei;
    JBBL jb1,jb2; 

    SHBL sh;

    //	OctCell *leafp;
    OctCell *pincell=NULL;
    vector< OctCell * > grid;


    ghostpt.x = ghostcell->xc1;
    ghostpt.y = ghostcell->yc1;	
    if(in2out==0){
        jb2=State_l;

        jbblToShbl(&jb2, &sh);

        ghostcell->dof[0][0] = sh.q;
        ghostcell->dof[1][0]= sh.qu;
        ghostcell->dof[2][0]= sh.qv;
        ghostcell->dof[3][0]= sh.te;
        for( i1=0; i1<4; ++i1) 
            for(j1=1; j1<nDOF; ++j1)
                ghostcell->dof[i1][j1]=0.0;

    }
    else {
        symtrcpt.x=2.0*xOutlet-ghostpt.x;
        symtrcpt.y=ghostpt.y;
        nx=-1.0;
        ny=0.0;
        getGridIndex2(ghostcell, &symtrcpt, pincell);
        findIntpltCell(pincell, ghostpt, symtrcpt, grid);
        chazhi4smtrictPshbl(grid, symtrcpt, jb1); 
        jbblToShbl(&jb1, &sh);

        ghostcell->dof[0][0] = sh.q;
        ghostcell->dof[1][0]= sh.qu;
        ghostcell->dof[2][0]= sh.qv;
        ghostcell->dof[3][0]= sh.te;
        for( i1=0; i1<4; ++i1) 
            for(j1=1; j1<nDOF; ++j1)
                ghostcell->dof[i1][j1]=0.0;
    }
}
extern vector<OctCell *> InletBndry;
extern vector<OctCell *> OutletBndry;
void inlet_outletBndry()
{
    for(vector<OctCell *> ::size_type isz=0; isz<InletBndry.size();++isz){
        getInletOutletBndry(InletBndry[isz], 0);
    }

    for(vector<OctCell *> ::size_type isz=0; isz<OutletBndry.size();++isz){
        getInletOutletBndry(OutletBndry[isz], 1);
    }
}

void farFieldBoundary()
{
    int j,k;
    int i1, j1;
    double nx;
    double ny;
    SHBL sh;
    JBBL jb1, jb2;

    double speed, a, vbn;

    vbn=0.0;
    //-----------------------------------------------------------------------

    //----远场边界条件，位于和y坐标轴垂直的下部面----begin--	 
    nx=0.0;
    ny = 1.0;

    for(k=1; k <= Nx; k++)
    {
        //----下部那个面-------begin--
        sh.q = bodygrid[2*Nx +k].dof[0][0];
        sh.qu = bodygrid[2*Nx +k].dof[1][0];
        sh.qv = bodygrid[2*Nx +k].dof[2][0];
        sh.te = bodygrid[2*Nx +k].dof[3][0];

        shblToJbbl(&sh, &jb1);

        a = sqrt(GAMMA * jb1.p / jb1.q);

        speed = jb1.u * nx + jb1.v * ny;

        if(speed > 0.0 && fabs(speed) >= a)  //inflow
        {
            jb2.q=initQ;
            jb2.p=initP;
            jb2.u=uFree;
            jb2.v=vFree;

            jbblToShbl(&jb2, &sh);

            bodygrid[Nx +k].dof[0][0] = sh.q;
            bodygrid[Nx +k].dof[1][0]= sh.qu;
            bodygrid[Nx +k].dof[2][0]= sh.qv;
            bodygrid[Nx +k].dof[3][0]= sh.te;
            for( i1=0; i1<4; ++i1) 
                for(j1=1; j1<nDOF; ++j1)
                    bodygrid[Nx +k].dof[i1][j1]=0.0;

        }
        else if(speed > 0.0 && fabs(speed) < a) //inflow
        {

            jb2.p = 0.5 * (initP + jb1.p) + 0.5 * jb1.q * a * (nx * (uFree - jb1.u) + ny * (vFree - jb1.v));
            vbn = 0.5 / (jb1.q * a) *  (initP - jb1.p) + 0.5 * (nx * (uFree + jb1.u) + ny * (vFree + jb1.v));

            jb2.q = initQ * pow(jb2.p/initP, 1./GAMMA);

            jb2.u = uFree + (vbn - (uFree * nx + vFree * ny )) * nx;
            jb2.v = vFree + (vbn - (uFree * nx + vFree * ny )) * ny;

            jbblToShbl(&jb2, &sh);

            bodygrid[Nx +k].dof[0][0] = sh.q;
            bodygrid[Nx +k].dof[1][0]= sh.qu;
            bodygrid[Nx +k].dof[2][0]= sh.qv;
            bodygrid[Nx +k].dof[3][0]= sh.te;
            for(i1=0; i1<4; ++i1) 
                for(j1=1; j1<nDOF; ++j1)
                    bodygrid[Nx +k].dof[i1][j1]=0.0;

        }
        else if(speed <= 0.0 && fabs(speed) >= a) //outflow
        {
            jb2.q=jb1.q;
            jb2.p=jb1.p;
            jb2.u=jb1.u;
            jb2.v=jb1.v;

            jbblToShbl(&jb2, &sh);

            bodygrid[Nx +k].dof[0][0] = sh.q;
            bodygrid[Nx +k].dof[1][0]= sh.qu;
            bodygrid[Nx +k].dof[2][0]= sh.qv;
            bodygrid[Nx +k].dof[3][0]= sh.te;
            for(i1=0; i1<4; ++i1) 
                for(j1=1; j1<nDOF; ++j1)
                    bodygrid[Nx +k].dof[i1][j1]=0.0;

        }
        else if(speed <= 0.0 && fabs(speed) < a) //outflow
        {

            jb2.p = 0.5 * (initP + jb1.p) + 0.5 * jb1.q * a * (nx * (uFree - jb1.u) + ny * (vFree - jb1.v));
            vbn = 0.5 / (jb1.q * a) *  (initP - jb1.p) + 0.5 * (nx * (uFree + jb1.u) + ny * (vFree + jb1.v));
            jb2.q = jb1.q * pow(jb2.p/jb1.p, 1./GAMMA);

            jb2.u = jb1.u + (vbn - (jb1.u * nx + jb1.v * ny)) * nx;

            jb2.v = jb1.v + (vbn - (jb1.u * nx + jb1.v * ny)) * ny;


            jbblToShbl(&jb2, &sh);

            bodygrid[Nx +k].dof[0][0] = sh.q;
            bodygrid[Nx +k].dof[1][0]= sh.qu;
            bodygrid[Nx +k].dof[2][0]= sh.qv;
            bodygrid[Nx +k].dof[3][0]= sh.te;
            for(i1=0; i1<4; ++i1) 
                for( j1=1; j1<nDOF; ++j1)
                    bodygrid[Nx +k].dof[i1][j1]=0.0;

        }

        //---虚拟网格外插------
        for(i1=0; i1<4; ++i1) {
            bodygrid[k].dof[i1][0] = 2.0*bodygrid[Nx+k].dof[i1][0]
                - bodygrid[2*Nx+k].dof[i1][0];
            for(j1=1; j1<nDOF;++j1)
                bodygrid[k].dof[i1][j1] = 0.0;
        }

        //----下部那个面-------end--
    }


    //////////////////////////////////////////////////////////////////////////////////---------------------2008.7.1

    //     for(j=1;j<=Ny;j++) {
    //       for(k=1;k<=Nx;k++) {	         				       		       
    //			bodygrid[(j-1)*Nx +k]

    //----远场边界条件，位于和y坐标轴垂直的上部面----begin--	 
    nx=0.0;
    ny = -1.0;


    for(k=1; k <= Nx; k++)
    {
        //----上部那个面-------begin--
        sh.q = bodygrid[(Ny-3)*Nx+k].dof[0][0];
        sh.qu = bodygrid[(Ny-3)*Nx+k].dof[1][0];
        sh.qv = bodygrid[(Ny-3)*Nx+k].dof[2][0];
        sh.te = bodygrid[(Ny-3)*Nx+k].dof[3][0];

        shblToJbbl(&sh, &jb1);

        a = sqrt(GAMMA * jb1.p / jb1.q);

        speed = jb1.u * nx + jb1.v * ny;

        if(speed > 0.0 && fabs(speed) >= a)  //inflow
        {
            jb2.q=initQ;
            jb2.p=initP;
            jb2.u=uFree;
            jb2.v=vFree;

            jbblToShbl(&jb2, &sh);

            bodygrid[(Ny-2)*Nx+k].dof[0][0] = sh.q;
            bodygrid[(Ny-2)*Nx+k].dof[1][0]= sh.qu;
            bodygrid[(Ny-2)*Nx+k].dof[2][0] = sh.qv;
            bodygrid[(Ny-2)*Nx+k].dof[3][0] = sh.te;

            for(i1=0; i1<4; ++i1)
                for(j1=1; j1<nDOF; ++j1) 
                    bodygrid[(Ny-2)*Nx+k].dof[i1][j1] = 0.0;

        }
        else if(speed > 0.0 && fabs(speed) < a)  //inflow
        {


            jb2.p = 0.5 * (initP + jb1.p) + 0.5 * jb1.q * a * (nx * (uFree - jb1.u) + ny * (vFree - jb1.v));
            vbn = 0.5 / (jb1.q * a) *  (initP - jb1.p) + 0.5 * (nx * (uFree + jb1.u) + ny * (vFree + jb1.v));
            jb2.q = initQ * pow(jb2.p/initP, 1./GAMMA);
            jb2.u = uFree + (vbn - (uFree * nx + vFree * ny )) * nx;
            jb2.v = vFree + (vbn - (uFree * nx + vFree * ny )) * ny;

            jbblToShbl(&jb2, &sh);

            bodygrid[(Ny-2)*Nx+k].dof[0][0] = sh.q;
            bodygrid[(Ny-2)*Nx+k].dof[1][0]= sh.qu;
            bodygrid[(Ny-2)*Nx+k].dof[2][0] = sh.qv;
            bodygrid[(Ny-2)*Nx+k].dof[3][0] = sh.te;

            for(i1=0; i1<4; ++i1)
                for(j1=1; j1<nDOF; ++j1) 
                    bodygrid[(Ny-2)*Nx+k].dof[i1][j1] = 0.0;
        }
        else if(speed <= 0.0 && fabs(speed) >= a)
        {
            jb2.q=jb1.q;
            jb2.p=jb1.p;
            jb2.u=jb1.u;
            jb2.v=jb1.v;

            jbblToShbl(&jb2, &sh);

            bodygrid[(Ny-2)*Nx+k].dof[0][0] = sh.q;
            bodygrid[(Ny-2)*Nx+k].dof[1][0]= sh.qu;
            bodygrid[(Ny-2)*Nx+k].dof[2][0] = sh.qv;
            bodygrid[(Ny-2)*Nx+k].dof[3][0] = sh.te;

            for(i1=0; i1<4; ++i1)
                for(j1=1; j1<nDOF; ++j1) 
                    bodygrid[(Ny-2)*Nx+k].dof[i1][j1] = 0.0;
        }
        else if(speed <= 0.0 && fabs(speed) < a)  //outflow
        {
            jb2.p = 0.5 * (initP + jb1.p) + 0.5 * jb1.q * a * (nx * (uFree - jb1.u) + ny * (vFree - jb1.v));
            vbn = 0.5 / (jb1.q * a) *  (initP - jb1.p) + 0.5 * (nx * (uFree + jb1.u) + ny * (vFree + jb1.v));
            jb2.q = jb1.q * pow(jb2.p/jb1.p, 1./GAMMA);

            jb2.u = jb1.u + (vbn - (jb1.u * nx + jb1.v * ny)) * nx;

            jb2.v = jb1.v + (vbn - (jb1.u * nx + jb1.v * ny)) * ny;

            jbblToShbl(&jb2, &sh);

            bodygrid[(Ny-2)*Nx+k].dof[0][0] = sh.q;
            bodygrid[(Ny-2)*Nx+k].dof[1][0]= sh.qu;
            bodygrid[(Ny-2)*Nx+k].dof[2][0] = sh.qv;
            bodygrid[(Ny-2)*Nx+k].dof[3][0] = sh.te;

            for(i1=0; i1<4; ++i1)
                for(j1=1; j1<nDOF; ++j1) 
                    bodygrid[(Ny-2)*Nx+k].dof[i1][j1] = 0.0;
        }

        //--虚拟网格外插------
        for (i1=0; i1<4; ++i1) {
            bodygrid[(Ny-1)*Nx+k].dof[i1][0] = 2.0*bodygrid[(Ny-2)*Nx+k].dof[i1][0]
                - bodygrid[(Ny-3)*Nx+k].dof[i1][0];
            for(j1=1; j1<nDOF; ++j1)
                bodygrid[(Ny-1)*Nx+k].dof[i1][j1] = 0.0;
        }
        //----------------------------------------------
        //----上部那个面-------end--
    }


    //-------------------右侧面---------0000000000000------begin---
    nx = -1.0;
    ny=0.0;

    for(j=1; j <= Ny; j++)
    {
        /*-----------begin--*/
        sh.q = bodygrid[(j-1)*Nx+Nx-2].dof[0][0];
        sh.qu = bodygrid[(j-1)*Nx+Nx-2].dof[1][0];
        sh.qv = bodygrid[(j-1)*Nx+Nx-2].dof[2][0];
        sh.te = bodygrid[(j-1)*Nx+Nx-2].dof[3][0];

        shblToJbbl(&sh, &jb1);

        a = sqrt(GAMMA * jb1.p / jb1.q);

        speed = jb1.u * nx + jb1.v * ny;

        if(speed > 0.0 && fabs(speed) >= a)
        {
            jb2.q=initQ;
            jb2.p=initP;
            jb2.u=uFree;
            jb2.v=vFree;

            jbblToShbl(&jb2, &sh);

            bodygrid[(j-1)*Nx+Nx-1].dof[0][0]= sh.q;
            bodygrid[(j-1)*Nx+Nx-1].dof[1][0] = sh.qu;
            bodygrid[(j-1)*Nx+Nx-1].dof[2][0] = sh.qv;
            bodygrid[(j-1)*Nx+Nx-1].dof[3][0] = sh.te;

            for(i1=0; i1<4; ++i1)
                for(j1=1; j1<nDOF; ++j1)
                    bodygrid[(j-1)*Nx+Nx-1].dof[i1][j1]=0.0;
        }
        else if(speed > 0.0 && fabs(speed) < a)  //inflow
        {

            jb2.p = 0.5 * (initP + jb1.p) + 0.5 * jb1.q * a * (nx * (uFree - jb1.u) + ny * (vFree - jb1.v));
            vbn = 0.5 / (jb1.q * a) * (initP - jb1.p) + 0.5 * (nx * (uFree + jb1.u) + ny * (vFree + jb1.v));

            jb2.q = initQ * pow(jb2.p/initP, 1./GAMMA);

            jb2.u = uFree + (vbn - (uFree * nx + vFree * ny)) * nx;
            jb2.v = vFree + (vbn - (uFree * nx + vFree * ny)) * ny;

            jbblToShbl(&jb2, &sh);

            bodygrid[(j-1)*Nx+Nx-1].dof[0][0]= sh.q;
            bodygrid[(j-1)*Nx+Nx-1].dof[1][0] = sh.qu;
            bodygrid[(j-1)*Nx+Nx-1].dof[2][0] = sh.qv;
            bodygrid[(j-1)*Nx+Nx-1].dof[3][0] = sh.te;

            for(i1=0; i1<4; ++i1)
                for(j1=1; j1<nDOF; ++j1)
                    bodygrid[(j-1)*Nx+Nx-1].dof[i1][j1]=0.0;
        }
        else if(speed <= 0.0 && fabs(speed) >= a)
        {
            jb2.q=jb1.q;
            jb2.p=jb1.p;
            jb2.u=jb1.u;
            jb2.v=jb1.v;

            jbblToShbl(&jb2, &sh);

            bodygrid[(j-1)*Nx+Nx-1].dof[0][0]= sh.q;
            bodygrid[(j-1)*Nx+Nx-1].dof[1][0] = sh.qu;
            bodygrid[(j-1)*Nx+Nx-1].dof[2][0] = sh.qv;
            bodygrid[(j-1)*Nx+Nx-1].dof[3][0] = sh.te;

            for(i1=0; i1<4; ++i1)
                for(j1=1; j1<nDOF; ++j1)
                    bodygrid[(j-1)*Nx+Nx-1].dof[i1][j1]=0.0;
        }
        else if(speed <= 0.0 && fabs(speed) < a)  //outflow
        {

            jb2.p = 0.5 * (initP + jb1.p) + 0.5 * jb1.q * a * (nx * (uFree - jb1.u) + ny * (vFree - jb1.v));
            vbn = 0.5 / (jb1.q * a) * (initP - jb1.p) + 0.5 * (nx * (uFree + jb1.u) + ny * (vFree + jb1.v));
            jb2.q = jb1.q * pow(jb2.p/jb1.p, 1./GAMMA);

            jb2.u = jb1.u + (vbn - (jb1.u * nx + jb1.v * ny)) * nx;

            jb2.v = jb1.v + (vbn - (jb1.u * nx + jb1.v * ny)) * ny;


            jbblToShbl(&jb2, &sh);

            bodygrid[(j-1)*Nx+Nx-1].dof[0][0]= sh.q;
            bodygrid[(j-1)*Nx+Nx-1].dof[1][0] = sh.qu;
            bodygrid[(j-1)*Nx+Nx-1].dof[2][0] = sh.qv;
            bodygrid[(j-1)*Nx+Nx-1].dof[3][0] = sh.te;

            for(i1=0; i1<4; ++i1)
                for(j1=1; j1<nDOF; ++j1)
                    bodygrid[(j-1)*Nx+Nx-1].dof[i1][j1]=0.0;
        }

        //---虚拟网格外插------
        for(i1=0; i1<4; ++i1) {
            bodygrid[(j-1)*Nx+Nx].dof[i1][0]
                = 2.0*bodygrid[(j-1)*Nx+Nx-1].dof[i1][0] - bodygrid[(j-1)*Nx+Nx-2].dof[i1][0];
            for(j1=1; j1<nDOF; ++j1)
                bodygrid[(j-1)*Nx+Nx].dof[i1][j1]=0.0;	
        }	       

    }

    //----右侧面------------------end--------

    //-------------------左侧面----begin------
    nx = 1.0;
    ny=0.0;

    for(j=1; j <= Ny; j++)
    {
        //-----------begin--------
        sh.q = bodygrid[(j-1)*Nx+3].dof[0][0];
        sh.qu = bodygrid[(j-1)*Nx+3].dof[1][0];
        sh.qv = bodygrid[(j-1)*Nx+3].dof[2][0];
        sh.te = bodygrid[(j-1)*Nx+3].dof[3][0];

        shblToJbbl(&sh, &jb1);

        a = sqrt(GAMMA * jb1.p / jb1.q);

        speed = jb1.u * nx + jb1.v * ny;

        if(speed > 0.0 && fabs(speed) >= a)
        {
            jb2.q=initQ;
            jb2.p=initP;
            jb2.u=uFree;
            jb2.v=vFree;

            jbblToShbl(&jb2, &sh);

            bodygrid[(j-1)*Nx+2].dof[0][0] = sh.q;
            bodygrid[(j-1)*Nx+2].dof[1][0]= sh.qu;
            bodygrid[(j-1)*Nx+2].dof[2][0] = sh.qv;
            bodygrid[(j-1)*Nx+2].dof[3][0]= sh.te;
            for(i1=0; i1<4; ++i1) 
                for(j1=1; j1<nDOF; ++j1)
                    bodygrid[(j-1)*Nx+2].dof[i1][j1]=0.0;

        }
        else if(speed > 0.0 && fabs(speed) < a)  //inflow
        {

            jb2.p = 0.5 * (initP + jb1.p) + 0.5 * jb1.q * a * (nx * (uFree - jb1.u) + ny * (vFree - jb1.v));
            vbn = 0.5 / (jb1.q * a) *  (initP - jb1.p) + 0.5 * (nx * (uFree + jb1.u) + ny * (vFree + jb1.v));

            jb2.q = initQ * pow(jb2.p/initP, 1./GAMMA);

            jb2.u = uFree + (vbn - (uFree * nx + vFree * ny)) * nx;
            jb2.v = vFree + (vbn - (uFree * nx + vFree * ny)) * ny;

            jbblToShbl(&jb2, &sh);

            bodygrid[(j-1)*Nx+2].dof[0][0] = sh.q;
            bodygrid[(j-1)*Nx+2].dof[1][0]= sh.qu;
            bodygrid[(j-1)*Nx+2].dof[2][0] = sh.qv;
            bodygrid[(j-1)*Nx+2].dof[3][0]= sh.te;
            for(i1=0; i1<4; ++i1) 
                for(j1=1; j1<nDOF; ++j1)
                    bodygrid[(j-1)*Nx+2].dof[i1][j1]=0.0;
        }
        else if(speed <= 0.0 && fabs(speed) >= a)
        {
            jb2.q=jb1.q;
            jb2.p=jb1.p;
            jb2.u=jb1.u;
            jb2.v=jb1.v;

            jbblToShbl(&jb2, &sh);

            bodygrid[(j-1)*Nx+2].dof[0][0] = sh.q;
            bodygrid[(j-1)*Nx+2].dof[1][0]= sh.qu;
            bodygrid[(j-1)*Nx+2].dof[2][0] = sh.qv;
            bodygrid[(j-1)*Nx+2].dof[3][0]= sh.te;
            for(i1=0; i1<4; ++i1) 
                for(j1=1; j1<nDOF; ++j1)
                    bodygrid[(j-1)*Nx+2].dof[i1][j1]=0.0;
        }
        else if(speed <= 0.0 && fabs(speed) < a)  //outflow
        {

            jb2.p = 0.5 * (initP + jb1.p) + 0.5 * jb1.q * a * (nx * (uFree - jb1.u) + ny * (vFree - jb1.v));
            vbn = 0.5 / (jb1.q * a) * (initP - jb1.p) + 0.5 * (nx * (uFree + jb1.u) + ny * (vFree + jb1.v));
            jb2.q = jb1.q * pow(jb2.p/jb1.p, 1./GAMMA);

            jb2.u = jb1.u + (vbn - (jb1.u * nx + jb1.v * ny)) * nx;

            jb2.v = jb1.v + (vbn - (jb1.u * nx + jb1.v * ny)) * ny;


            jbblToShbl(&jb2, &sh);

            bodygrid[(j-1)*Nx+2].dof[0][0] = sh.q;
            bodygrid[(j-1)*Nx+2].dof[1][0]= sh.qu;
            bodygrid[(j-1)*Nx+2].dof[2][0] = sh.qv;
            bodygrid[(j-1)*Nx+2].dof[3][0]= sh.te;
            for(i1=0; i1<4; ++i1) 
                for(j1=1; j1<nDOF; ++j1)
                    bodygrid[(j-1)*Nx+2].dof[i1][j1]=0.0;
        }

        //---虚拟网格外插------
        for(i1=0; i1<4; ++i1) {
            bodygrid[(j-1)*Nx+1].dof[i1][0]
                = 2.0*bodygrid[(j-1)*Nx+2].dof[i1][0] - bodygrid[(j-1)*Nx+3].dof[i1][0];
            for(j1=1; j1<nDOF; ++j1)
                bodygrid[(j-1)*Nx+1].dof[i1][j1]=0.0;
        }
    }

    //----左侧面-------------------------end----
}


/*
   void copy_jbbl(JBBL *jbnp, JBBL *jbtmpr);
   void copy_shbl(SHBL *jbnp, SHBL *jbtmpr);
   void farFieldBoundary0()
   {
   int iy, ix;

   JBBL jb2;
   SHBL sh2;


   jb2.q=initQ;
   jb2.p=initP;
   jb2.u=uFree;
   jb2.v=vFree;
   jb2.T=jb2.p/jb2.q /Rstar;
   jb2.k=initK;
   jb2.omg=initOmg;

   jbblToShbl(&jb2, &sh2);

   for(iy=0; iy<=Ny+1; iy++)
   {	
   for(ix=0; ix<=Nx+1; ix++)
   {
   if(iy == 0 || iy == 1 || iy == 2 
   || iy == Ny+1 || iy == Ny || iy == Ny-1 
   || ix == 0 || ix == 1 || ix == 2 
   || ix == Nx+1 || ix == Nx || ix == Nx-1 )
   {
   copy_jbbl(&jb2, &bodygrid[iy][ix].jbnp);
   copy_shbl(&sh2, &bodygrid[iy][ix].shbl1np);
   }
   }
   }
   }

*/
