//nihao 
//OctCell *NorthNeighbor(OctCell *pp);
//OctCell *WestNeighbor(OctCell *pp);
//OctCell *EastNeighbor(OctCell *pp);
//OctCell *SouthNeighbor(OctCell *pp);
#include"non_uniform_grid.h"
#include"findneighbor.h"
//using vim to edit the program 

extern OctCell *bodygrid;

double minmod3n(double a1, double a2, double a3)
{

    double tmp;
    if(a1>0.0 && a2>0.0 && a3>0.0) 
        tmp=min3(a1,a2,a3);
    else if(a1<0.0 && a2 < 0.0 && a3<0.0)
        tmp=max3(a1,a2,a3);
    else 
        tmp=0.0;
    return tmp;

}


double minmod3n_general(double a1, double a2, double a3, double dxy)
{
    double rt;
    if(fabs(a1) <= M4TVBM*dxy*dxy) rt = a1;
    else 
        rt=minmod3n(a1,a2,a3);
    return rt;
}


void sonCellAvg(int nson,OctCell *pcell, double dof0[4]);

void tvbm_limiter(OctCell *pcell)
{
    double dx, dy;
    double tmpx, tmpy;

    short int lev;

    double Rx[4][4]={0.0}; //Rx右特征矩阵
    double invRx[4][4]={0.0};  //左特征向量矩阵

    double Ry[4][4]={0.0}; //Ry右特征矩阵
    double invRy[4][4]={0.0};  

    dx=pcell->dx;
    dy=pcell->dy;

    lev=pcell->level;
    double dof[4][6]={0.0};

    double upxdof0[4]={0.0};
    double umxdof0[4]={0.0};
    double umydof0[4]={0.0};
    double upydof0[4]={0.0};

    OctCell *pcelle=NULL;
    OctCell *pcells=NULL;
    OctCell *pcelln=NULL;
    OctCell *pcellw=NULL;

    pcelle=EastNeighbor(pcell);
    pcells=SouthNeighbor(pcell);
    pcellw=WestNeighbor(pcell);
    pcelln=NorthNeighbor(pcell);

    int nson0,nson;//the number of the son in his parent

    for(int i=0;i<4;i++) 
        for(int j=0; j< nDOF; j++) dof[i][j]= pcell->dof[i][j];

    nson0=pcell->NOchildren;//该网格单元在子网格中的编号
    if(lev==pcelle->level){
        for(int i=0;i<4;++i){
            upxdof0[i]=pcelle->dof[i][0];
        }
    }
    else {
        switch(nson0) {
        case 1:
            nson=0;
            break;
        case 2:
            nson=3;
            break;
        }
        sonCellAvg(nson,pcelle,upxdof0);
    }

    if(lev==pcellw->level){
        for(int i=0;i<4;++i){
            umxdof0[i]=pcellw->dof[i][0];
        }
    }
    else {
        switch(nson0) {
        case 0:
            nson=1;
            break;
        case 3:
            nson=2;
            break;
        }
        sonCellAvg(nson,pcellw,umxdof0);
    }

    if(lev==pcells->level){
        for(int i=0;i<4;++i){
            umydof0[i]=pcells->dof[i][0];
        }
    }
    else {
        switch(nson0) {
        case 0:
            nson=3;
            break;
        case 1:
            nson=2;
            break;
        }
        sonCellAvg(nson,pcells,umydof0);
    }

    if(lev==pcelln->level){
        for(int i=0;i<4;++i){
            upydof0[i]=pcelln->dof[i][0];
        }
    }
    else {
        switch(nson0) {
        case 3:
            nson=0;
            break;
        case 2:
            nson=1;
            break;
        }
        sonCellAvg(nson,pcelln,upydof0);
    }

#ifdef TZBL_LIM

    double u, v, pre, a0, gma, han, a2;

    gma=GAMMA-1.0;

    u=dof[1][0]/dof[0][0];
    v=dof[2][0]/dof[0][0];

    pre=gma*(dof[3][0]-0.5*dof[0][0]*(u*u+v*v));

    a0=sqrt(GAMMA*pre/dof[0][0]);

    a2=a0*a0;

    han=(dof[3][0]+pre)/dof[0][0];

    //--qiu Rx RY, invRx, invRy
    /*
       double rcm, b1, b2, t0x, t1, t2, t3,qm, t0y, t0;

       qm=0.5*(u*u+v*v);
       rcm = 1.0 / a0;
       b1 = gma * rcm * rcm;
       b2 = qm * b1;
       t0x = u * rcm;
       t1 = b1 * u;
       t2 = 0.5 * b1;
       t3 = b1 * v;

       invRx[0][0]=0.5 * ( b2 + t0x );

       invRx[0][1]=-0.5 * ( t1 + rcm );
       invRx[0][2]=  -0.5 * t3;
       invRx[0][3]= t2;


       invRx[1][0]=  - v;    
       invRx[1][1]=  0.0;
       invRx[1][2]=  1.0;
       invRx[1][3]=  0.0;


       invRx[2][0]=  1.0 - b2;
       invRx[2][1]= t1;
       invRx[2][2]=t3;
       invRx[2][3]= -b1;

       invRx[3][0]=   0.5 * ( b2 - t0x );
       invRx[3][1]= -0.5 * ( t1 - rcm );
       invRx[3][2]= -0.5 *   t3;
       invRx[3][3]= t2;     

       t0y = v * rcm;
       t1 = b1 * v;
       t2 = 0.5 * b1;
       t3 = b1 * u;

       invRy[0][0] = 0.5 * ( b2 + t0y );
       invRy[0][1] = -0.5 * t3;
       invRy[0][2] = -0.5 * ( t1 + rcm );
       invRy[0][3] = t2;

       invRy[1][0] = - u;
       invRy[1][1] = 1.0;
       invRy[1][2] = 0.0;
       invRy[1][3] = 0.0;

       invRy[2][0] = 1. - b2;
       invRy[2][1] = t3;
       invRy[2][2] = t1;
       invRy[2][3]= -b1;

       invRy[3][0]=  0.5 * ( b2 - t0y );
       invRy[3][1] = -0.5 *   t3;
       invRy[3][2] = -0.5 * ( t1 - rcm );
       invRy[3][3] = t2 ;

       t0 = u * a0;

       Rx[0][0] = 1.0;
       Rx[0][1]  = 0.0;
       Rx[0][2]  = 1.0;
       Rx[0][3]  = 1.0;

       Rx[1][0] = u- a0;
       Rx[1][1] = 0.0;
       Rx[1][2] = u;
       Rx[1][3] = u + a0;

    Rx[2][0] = v;
    Rx[2][1]  = 1.0;
    Rx[2][2] = v;
    Rx[2][3]  = v;

    Rx[3][0]  = han - t0;
    Rx[3][1] = v;
    Rx[3][2] = qm;
    Rx[3][3] = han+ t0;

    t0 = v* a0;

    Ry[0][0] = 1.0;
    Ry[0][1] = 0.0;
    Ry[0][2] = 1.0;
    Ry[0][3] = 1.0;

    Ry[1][0] = u;
    Ry[1][1]  = 1.0;
    Ry[1][2]  = u;
    Ry[1][3]  = u;

    Ry[2][0]  = v - a0;
    Ry[2][1] = 0.0;
    Ry[2][2] = v;
    Ry[2][3] = v + a0;

    Ry[3][0] = han- t0;
    Ry[3][1] = u;
    Ry[3][2] = qm;
    Ry[3][3] = han + t0;
    */           

    //--qiu
    ///liu
        Rx[0][0]=1.0;
    Rx[0][1]=1.0;
    Rx[0][2]=0.0;
    Rx[0][3]=1.0;

    Rx[1][0]= u-a0;
    Rx[1][1]= u;
    Rx[1][2]= 0.0;
    Rx[1][3]= u+a0;

    Rx[2][0]= v;
    Rx[2][1]= v;
    Rx[2][2]= 1.0;
    Rx[2][3]= v;

    Rx[3][0]= han-u*a0;
    Rx[3][1]= 0.5*(u*u+v*v);
    Rx[3][2]= v;
    Rx[3][3]= han+u*a0;

    Ry[0][0]=1.0;
    Ry[0][1]=0.0;
    Ry[0][2]=1.0;
    Ry[0][3]=1.0;

    Ry[1][0]= u;
    Ry[1][1]= 1.0;
    Ry[1][2]= u;
    Ry[1][3]= u;

    Ry[2][0]= v-a0;
    Ry[2][1]= 0.0;
    Ry[2][2]= v;
    Ry[2][3]= v+a0;

    Ry[3][0]= han-v*a0;
    Ry[3][1]= u;
    Ry[3][2]= 0.5*(u*u+v*v);
    Ry[3][3]= han+v*a0;    

    invRx[0][0]=han+a0/gma*(u-a0);
    invRx[0][1]= -(u+a0/gma);
    invRx[0][2]=  -v;
    invRx[0][3]= 1.0;

    invRx[1][0]= -2.0*han+4.0*a2/gma;    
    invRx[1][1]= 2.0*u;
    invRx[1][2]= 2.0*v;
    invRx[1][3]= -2.0;

    invRx[2][0]= -2.0*v*a2/gma;
    invRx[2][1]= 0.0;
    invRx[2][2]=2.0*a2/gma;
    invRx[2][3]= 0.0;

    invRx[3][0]= han- a0*(u+a0)/gma;
    invRx[3][1]= -u+a0/gma;
    invRx[3][2]= -v;
    invRx[3][3]=1.0;

    invRy[0][0]= han+ a0/gma*(v-a0);
    invRy[0][1]= -u;
    invRy[0][2]= -(v+a0/gma);
    invRy[0][3]= 1.0;

    invRy[1][0]= -2.0*u*a2/gma;
    invRy[1][1]= 2.0*a2/gma;
    invRy[1][2]= 0.0;
    invRy[1][3]= 0.0;

    invRy[2][0]= -2.0*han+4.0*a2/gma; 
    invRy[2][1]= 2.0*u;
    invRy[2][2]= 2.0*v;
    invRy[2][3]= -2.0;

    invRy[3][0]= han- a0/gma*(v+a0);
    invRy[3][1]= -u;
    invRy[3][2]= -v+a0/gma;
    invRy[3][3]= 1.0;

    for(int i=0; i<4; i++) {
        for(int j=0; j<4; j++) {
            invRx[i][j] *= 0.5*gma/a2;
            invRy[i][j] *= 0.5*gma/a2;   
        }
    }   
    ///liu

    // double udx[4]={0.0}, udy[4]={0.0}; 
    double udxlim[4]={0.0}, udylim[4]={0.0};    

    double upx[4]={0.0}, umx[4]={0.0}, upy[4]={0.0}, umy[4]={0.0};
    //----------------------------------------------limt1

    double  ux1tm[4]={0.0}, ux2tm[4]={0.0},  uy1tm[4]={0.0}, uy2tm[4]={0.0};

    double ux1t[4]={0.0}, ux2t[4]={0.0},  uy1t[4]={0.0}, uy2t[4]={0.0};

    for(int i=0; i<4; i++) {
        for(int j=0; j<4; j++) {

            //           udx[i]+= invRx[i][j]*dof[j][1];

            //        udy[i]+= invRy[i][j]*dof[j][2];  

            ux1t[i]+= invRx[i][j]*(dof[j][1]+ 2.0/3.0*dof[j][4]- 1.0/3.0*dof[j][5]);
            ux2t[i]+= -invRx[i][j]*(-dof[j][1]+ 2.0/3.0*dof[j][4]- 1.0/3.0*dof[j][5]);           

            uy1t[i]+= invRy[i][j]*(dof[j][2]- 1.0/3.0*dof[j][4] + 2.0/3.0*dof[j][5]);
            uy2t[i]+= -invRy[i][j]*(-dof[j][2]- 1.0/3.0*dof[j][4]+ 2.0/3.0*dof[j][5]);   


            upx[i]+= invRx[i][j]*(upxdof0[j]-dof[j][0]);

            umx[i]+= invRx[i][j]*(dof[j][0]-umxdof0[j]);

            upy[i]+= invRy[i][j]*(upydof0[j]-dof[j][0]);
            umy[i]+= invRy[i][j]*(dof[j][0]-umydof0[j]);
        }
    }

    int nlimflg=0;

    for(int i=0; i<4; i++)  {     

        ux1tm[i]=minmod3n_general(ux1t[i],  upx[i],  umx[i], dx);

        ux2tm[i]=minmod3n_general(ux2t[i],  upx[i], umx[i], dx);

        uy1tm[i]=minmod3n_general(uy1t[i],  upy[i],  umy[i], dy);

        uy2tm[i]=minmod3n_general(uy2t[i],  upy[i], umy[i], dy);             


        if(fabs(ux1t[i]-ux1tm[i]) > ERRS  || fabs(ux2t[i]-ux2tm[i]) > ERRS ||
            fabs(uy1t[i]-uy1tm[i])> ERRS  || fabs(uy2t[i]-uy2tm[i])> ERRS ) ++nlimflg;
    }    


    double  rux1tm(0.0), rux2tm(0.0),  ruy1tm(0.0), ruy2tm(0.0);
    if(nlimflg!=0) {     
        for(int i=0; i<4; i++)  {
            rux1tm=rux2tm=ruy1tm=ruy2tm=0.0;

            for(int j=0; j<4; j++) {
                rux1tm+= Rx[i][j]*ux1tm[j];
                rux2tm+= Rx[i][j]*ux2tm[j];
                ruy1tm+= Ry[i][j]*uy1tm[j];
                ruy2tm+= Ry[i][j]*uy2tm[j];                
            }

            pcell->dof[i][1]=0.5*(rux1tm+rux2tm);
            pcell->dof[i][2]=0.5*(ruy1tm+ruy2tm);          
            pcell->dof[i][4]=rux1tm-rux2tm+0.5*(ruy1tm-ruy2tm);

            pcell->dof[i][5]=0.5*(rux1tm-rux2tm)+ruy1tm-ruy2tm;          

            pcell->dof[i][3]=0.0;
        } 

    }      

    //----------------------------------------lmit1     

    //limt2
    /*
       for(i=0; i<4; i++) {
       for(j=0; j<4; j++) {
       udx[i]+= invRx[i][j]*dof[j][1];            

       udy[i]+= invRy[i][j]*dof[j][2];

       upx[i]+= invRx[i][j]*(bodygrid[iy][ix+1].dof[j][0]-dof[j][0]);

       umx[i]+= invRx[i][j]*(dof[j][0]-bodygrid[iy][ix-1].dof[j][0]);

       upy[i]+= invRy[i][j]*(bodygrid[iy+1][ix].dof[j][0]-dof[j][0]);
       umy[i]+= invRy[i][j]*(dof[j][0]-bodygrid[iy-1][ix].dof[j][0]);
       }
       }

       for(i=0; i<4; i++)  {
    //对于方程组需要对特征变量使用限制器 

    tmpx= minmod3n_general(udx[i],  lim_muscl*upx[i],  lim_muscl*umx[i], dx);

    tmpy= minmod3n_general(udy[i],  lim_muscl*upy[i],  lim_muscl*umy[i], dy);

    if(fabs(udx[i]-tmpx)<= ERRS  && fabs(udy[i]-tmpy)<= ERRS)
    {	   
    bodygrid[iy][ix].lmtflg[i]=0;
    udxlim[i]=udx[i];
    udylim[i]=udy[i];
    }
    else { 
    bodygrid[iy][ix].lmtflg[i]=1;

    udxlim[i]=tmpx;
    udylim[i]=tmpy;

    bodygrid[iy][ix].dof[i][3] = bodygrid[iy][ix].dof[i][4]= bodygrid[iy][ix].dof[i][5]=0.0;
    //     bodygrid[iy][ix].dof0[i][3] = bodygrid[iy][ix].dof0[i][4]= bodygrid[iy][ix].dof0[i][5]=0.0;         
    }
    }  


    //-------------------------------------------------  limt2

    double rux[4]={0.0}, ruy[4]={0.0};

    for(i=0; i<4; i++) {
    if(bodygrid[iy][ix].lmtflg[i]==1) {
    for(j=0; j<4; j++) {
    rux[i]+= Rx[i][j]*udxlim[j];
    ruy[i]+= Ry[i][j]*udylim[j];
    }

    bodygrid[iy][ix].dof[i][1]=rux[i];
    bodygrid[iy][ix].dof[i][2]=ruy[i];

    //    bodygrid[iy][ix].dof[i][3]= bodygrid[iy][ix].dof[i][4]= bodygrid[iy][ix].dof[i][5]=0.0;            
    }
    }    

*/
#endif   

}


extern Node *HeadListAllGrid;
//边界网格使用限制器有问题，这里边界外的单元都设为0，会有影响 TODO
void tvdm_limiter4all()
{
    Node *current;
    current = HeadListAllGrid;  
    OctCell *pcell0;
    while(current != NULL)
    {		
        pcell0 = current->cell;
        if(pcell0->flag == 0 && current->flg<=2 ) {
#ifdef KXRCF_LIM
            if(pcell0->trb==true)
#endif
            tvbm_limiter(pcell0);
        }        
        current = current->next;
    }
}

//直接x,y是原始计算区域的坐标 [x_{i-1/2} x_{i+1/2}] 
double n_approx_sol(double x,double y, double xc1, double yc1, double dx, double dy, double *doftmp, int n );
extern double epsPos;
void pos_limiter(OctCell *pc)
{
    double dx,dy,xc,yc;
    double dof[4][6]={0.0};
    double cgsx[3];
    double cgsy[3];
    
    double cgslbx[4]={0.0};
    double cgslby[4]={0.0};

    double gslbwei[4]={0.0};
    double gslbpt[4]={0.0};
    double qmin=0.0, rtp,tht1,tht2, pre0;
    double u, v, gma1;
    double pre1;
    double ush[4],ujb[4];

    dx=pc->dx;
    dy=pc->dy;
    xc=pc->xc1;
    yc=pc->yc1;

    for(int i=0;i<4;i++) 
        for(int j=0; j< nDOF; j++) dof[i][j]=pc->dof[i][j];


    if(NLBpt==3){
        for(int i=0; i<3;++i){
            gslbwei[i]=weiLB3[i];
            gslbpt[i]=gsLB3pt[i];
        }
    }
    else{ 
        for(int i=0; i<4;++i){
            gslbwei[i]=weiLB4[i];
            gslbpt[i]=gsLB4pt[i];
        }
    }

    for(int i=0; i<3;++i){
        cgsx[i]=xc+0.5*dx*gspt[i];
        cgsy[i]=yc+0.5*dy*gspt[i];
    }

    for(int i=0; i<NLBpt;++i){
        cgslbx[i]=xc+0.5*dx*gslbpt[i];
        cgslby[i]=yc+0.5*dy*gslbpt[i];
    }

//limit density
    if(dof[0][0]<epsPos){
    //    dof[0][0]=epsPos;
        for(int mm=1;mm<nDOF;++mm)
            dof[0][mm]=0.0;
    }
    else{
        qmin=1.0e20;
        for(int i=0;i<3;++i){
            for(int j=0;j<NLBpt;++j){
                rtp=n_approx_sol(cgslbx[j],cgsy[i],
                    xc,yc,dx,dy, dof[0],nDOF);
                if(rtp<qmin) qmin=rtp;
            }
        }

        for(int i=0;i<3;++i){
            for(int j=0;j<NLBpt;++j){
                rtp=n_approx_sol(cgsx[i],cgslby[j],
                    xc,yc,dx,dy, dof[0],nDOF);
                if(rtp<qmin) qmin=rtp;
            }
        }

        tht1=min2(fabs((dof[0][0]-epsPos)/(dof[0][0]-qmin)), 1.0);

        for(int i=1;i<nDOF;++i)
            dof[0][i]*=tht1;
    }


    gma1=GAMMA-1.0;

    u=dof[1][0]/dof[0][0];
    v=dof[2][0]/dof[0][0];

    pre0=gma1*(dof[3][0]-0.5*dof[0][0]*(u*u+v*v));

    if(pre0<epsPos){
        for(int i=0;i<4;++i){
            for(int j=1;j<nDOF;j++){
                dof[i][j]=0.0;
            }
        }
    }
    else{
        tht2=1.0e20;
        for(int i=0;i<3;++i){
            for(int j=0;j<NLBpt;++j){
                for(int a=0;a<4;++a){
                    ush[a]=n_approx_sol(cgsx[i],cgslby[j],
                        xc,yc,dx,dy, dof[a],nDOF);
                }
                shbl2jbbl(ush,ujb);
                pre1=ujb[3];
                if(pre1>epsPos){
                    rtp=1.0;
                }
                else{
                    rtp=(pre0-epsPos)/(pre0-pre1);
                }

                if(rtp<tht2) tht2=rtp;
            }
        }


        for(int i=0;i<3;++i){
            for(int j=0;j<NLBpt;++j){
                for(int a=0;a<4;++a){
                    ush[a]=n_approx_sol(cgslbx[j],cgsy[i],
                        xc,yc,dx,dy, dof[a],nDOF);
                }
                shbl2jbbl(ush,ujb);
                pre1=ujb[3];
                if(pre1>epsPos){
                    rtp=1.0;
                }
                else{
                    rtp=(pre0-epsPos)/(pre0-pre1);
                }

                if(rtp<tht2) tht2=rtp;
            }
        }

        for(int j=0;j<4;++j)
            for(int i=1;i<nDOF;++i)
                dof[j][i]*=tht2;
    }

    for(int i=0;i<4;i++) 
        for(int j=1; j< nDOF; j++) pc->dof[i][j]= dof[i][j];
}

    
void pos_limiter4all()
{
    Node *current;
    current = HeadListAllGrid;  
    OctCell *pcell0;
    while(current != NULL)
    {		
        pcell0 = current->cell;
        if(pcell0->flag == 0 && current->flg<=2 ) {
            pos_limiter(pcell0);
        }        
        current = current->next;
    }
}















