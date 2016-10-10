
#include"non_uniform_grid.h"


//extern OctCell *bodygrid;
extern Node *HeadListAllGrid;
extern double dt;
extern double timesum;

void shblToJbbl(SHBL *, JBBL *);
void jbblToShbl(JBBL *, SHBL *);
void spatial_discretization_all();
void tvdm_limiter4all();
void face_flux_all();
void cal_upcx_upcy();

void inlet_outletBndry();
void farFieldBoundary();
void farFieldBoundary(int);
void wall_treatment(); 
void high_order_wall_treatment();
void tree_avg();
void tree_0bAvg();

void setAllLeafbAvgFalse();
const int NpreTimeStp=-25;
const double small_dt=5.0e-4;
//放到Runge-kutta步中   
//
/*************************************************************
  第一步RK之前先调用dof2dof0函数 leave the previous value.
 **************************************************************/
//void dof2dof0()  {
//   int ix,iy;
//   int imm, jmm;

//   for(iy=0; iy<= Ny+1; iy++) {
//       for(ix=0; ix <= Nx+1; ix++) {
//            for(imm=0; imm<4; imm++) 
//                for(jmm=0; jmm<nDOF; jmm++)
//                    bodygrid[iy][ix].dof0[imm][jmm]=bodygrid[iy][ix].dof[imm][jmm];
//        }
//    }
//}

void dof2dof0()
{
    Node *current;
    current = HeadListAllGrid;  
    OctCell *pcell0;
    while(current != NULL)
    {		
        pcell0 = current->cell;
        if(pcell0->flag % 2 == 0 ) {
            for(int im=0;im<4; ++im)
                for(int jm=0;jm<nDOF; ++jm)
                    pcell0->dof0[im][jm]=pcell0->dof[im][jm];
        }        
        current = current->next;
    }
}  			

void fn_RKTVD(OctCell *unp, int nstp, int rkstp)
{
    int im,jm;
#ifdef LOCAL_TIME 
    double dt=unp->dt1;
#endif
    if(nstp<NpreTimeStp){
        dt=small_dt;
    }

    switch(rkstp){
    case 0:
        for(im=0; im<4; im++) {
            for(jm=0; jm<nDOF; jm++) {
                unp->dof[im][jm]=unp->dof0[im][jm]
                    - dt*unp->invM[jm] * unp->residual[im][jm];
            }
        }
        break;
    case 1:
        for(im=0; im<4; im++) {
            for(jm=0; jm<nDOF; jm++) {
                unp->dof[im][jm]=0.75*unp->dof0[im][jm]
                    + 0.25*(unp->dof[im][jm]- dt*unp->invM[jm]
                        * unp->residual[im][jm]);
            }
        }
        break;
    case 2:
        for(im=0; im<4; im++) {
            for(jm=0; jm<nDOF; jm++) {
                unp->dof[im][jm]=1.0/3.0*unp->dof0[im][jm] 
                    + 2.0/3.0*(unp->dof[im][jm]
                        - dt*unp->invM[jm] * unp->residual[im][jm]);
            }
        }
        break;
    }
    unp->bAvg=true;
}

void fn_RKTVDAll(int tstp, int rkstp)
{
    Node *current;
    current = HeadListAllGrid;  
    OctCell *pcell0;
    while(current != NULL)
    {		
        pcell0 = current->cell;
        if(pcell0->flag == 0 && current->flg<=2 ) {
            fn_RKTVD(pcell0,tstp, rkstp);
        }        
        current = current->next;
    }
}

void setExtWallBndry();
void resd_smoothing();
void kxrcf_relate();
void set_kxrcf();
void face_kxrcf_all( );
void pos_limiter4all();
void runge_kutta(int tstp)
{	
    //step1	 
    dof2dof0();         		
    for(int rkstp=0;rkstp<3;++rkstp){
        if(Limiter01 != 0){
#ifdef KXRCF_LIM
            kxrcf_relate();
            face_kxrcf_all();
            set_kxrcf();
#endif
            tree_avg();
            tvdm_limiter4all();  

            setAllLeafbAvgFalse();
            tree_0bAvg();
        }
#if PosLimiter==1
        pos_limiter4all();
#endif

#if IS_TUBE==0
        farFieldBoundary();
#endif
        wall_treatment(); //放到Runge-kutta步中
#if IS_TUBE ==1
        setExtWallBndry();
        inlet_outletBndry();
#endif
        face_flux_all();

        spatial_discretization_all();
        fn_RKTVDAll(tstp,rkstp);

    }
}

/*
void RKTVDfirst(OctCell *unp, int nstp)
{
    int im,jm;
#ifdef LOCAL_TIME 
    double dt=unp->dt1;
#endif
    if(nstp<NpreTimeStp){
        dt=small_dt;
    }
    for(im=0; im<4; im++) {
        for(jm=0; jm<nDOF; jm++) {
            unp->dof[im][jm]=unp->dof0[im][jm]
                - dt*unp->invM[jm] * unp->residual[im][jm];
        }
    }
    unp->bAvg=true;
}


void RKTVDsecond(OctCell *unp,int nstp)
{
    int im,jm;
#ifdef LOCAL_TIME 
    double dt=unp->dt1;
#endif
    if(nstp<NpreTimeStp){
        dt=small_dt;
    }
    for(im=0; im<4; im++) {
        for(jm=0; jm<nDOF; jm++) {
            unp->dof[im][jm]=0.75*unp->dof0[im][jm]
                + 0.25*(unp->dof[im][jm]- dt*unp->invM[jm]
                    * unp->residual[im][jm]);
        }
    }
    unp->bAvg=true;
}


void RKTVDthird(OctCell *unp,int nstp)
{
    int im,jm;
#ifdef LOCAL_TIME 
    double dt=unp->dt1;
#endif
    if(nstp<NpreTimeStp){
        dt=small_dt;
    }
    for(im=0; im<4; im++) {
        for(jm=0; jm<nDOF; jm++) {
            unp->dof[im][jm]=1.0/3.0*unp->dof0[im][jm] 
                + 2.0/3.0*(unp->dof[im][jm]
                    - dt*unp->invM[jm] * unp->residual[im][jm]);
        }
    }
    unp->bAvg=true;
}

void RKTVDfirstAll(int tstp)
{
    Node *current;
    current = HeadListAllGrid;  
    OctCell *pcell0;
    while(current != NULL)
    {		
        pcell0 = current->cell;
        if(pcell0->flag == 0 && current->flg<=2 ) {
            RKTVDfirst(pcell0,tstp);
        }        
        current = current->next;
    }

}


void RKTVDsecondAll(int tstp)
{
    Node *current;
    current = HeadListAllGrid;  
    OctCell *pcell0;
    while(current != NULL)
    {		
        pcell0 = current->cell;
        if(pcell0->flag == 0 && current->flg<=2 ) {
            RKTVDsecond(pcell0,tstp);
        }        
        current = current->next;
    }

}

void RKTVDthirdAll(int tstp)
{
    Node *current;
    current = HeadListAllGrid;  
    OctCell *pcell0;
    while(current != NULL)
    {		
        pcell0 = current->cell;
        if(pcell0->flag == 0 && current->flg<=2 ) {
            RKTVDthird(pcell0,tstp);
        }        
        current = current->next;
    }

}

void runge_kutta(int tstp)
{	
    //step1	 
    dof2dof0();         		

    farFieldBoundary();

    //    cout<<"farfieldboundary"<<endl;
    //   #ifdef DUICHEN  
    //  cout<<"==="<<tstp<<endl;
    wall_treatment(); //放到Runge-kutta步中
    // cout<<"+++++"<<endl;
    //   #elif defined ILWILW
    //	 high_order_wall_treatment();
    //  #endif
    //	cout<<"wall treatment"<<endl;
    //

    //#ifdef LLFFLUX
    //       cal_upcx_upcy();
    //#endif

    face_flux_all();
    //   cout<<"hello facfluc\n";
    spatial_discretization_all();
    RKTVDfirstAll(tstp);
    if(Limiter01 != 0){
        tree_avg();
        tvdm_limiter4all();  

        setAllLeafbAvgFalse();
        tree_0bAvg();
    }
    //step2            
    farFieldBoundary();
    //  #ifdef DUICHEN  
    wall_treatment(); //放到Runge-kutta步中
    //  #elif defined ILWILW
    //	 high_order_wall_treatment();
    //  #endif   


    //#ifdef LLFFLUX
    //        cal_upcx_upcy();
    //#endif       
    face_flux_all();
    spatial_discretization_all();
    RKTVDsecondAll(tstp);	
    if(Limiter01 != 0) {
        tree_avg();
        tvdm_limiter4all();  

        setAllLeafbAvgFalse();
        tree_0bAvg();
    }
    //step3           
    farFieldBoundary();
    //   #ifdef DUICHEN  
    wall_treatment(); //放到Runge-kutta步中
    //   #elif defined ILWILW
    //	 high_order_wall_treatment();
    //   #endif

    //#ifdef LLFFLUX
    //        cal_upcx_upcy();
    //#endif
    face_flux_all();
    spatial_discretization_all();
    RKTVDthirdAll(tstp);  	
    if(Limiter01 != 0) {
        tree_avg();
        tvdm_limiter4all();

        setAllLeafbAvgFalse();
        tree_0bAvg();
    }
}*/




/*-------------------------------------------------------------------
  Enthalpy damping
  extern REAL uFree,vFree,initP,initQ;
  void enthalpyDamping()
  {
  REAL Cp, tH,u,v, tHinf;
  int iy,ix;
  REAL sgm=1.0e-3;
  Cp=GAMMA/(GAMMA-1.0) * Rstar;
  SHBL shu;
  JBBL jbu;

  REAL pls;
  tHinf = initT * Cp + 0.5*(uFree*uFree+vFree*vFree);

  for(iy=3; iy<=Ny-2; iy++)
  {
  for(ix=3; ix <= Nx-2; ix++)
  {
  if(bodygrid[iy][ix].flag == 0)
  {
  u=bodygrid[iy][ix].jbnp.u;
  v=bodygrid[iy][ix].jbnp.v;
  tH=Cp*bodygrid[iy][ix].jbnp.T + 0.5*(u*u+v*v);

  pls=1.0/(1.0+sgm*(tH-tHinf));

  bodygrid[iy][ix].shbl1np.q=pls*bodygrid[iy][ix].shbl1np.q;
  bodygrid[iy][ix].shbl1np.qu=pls*bodygrid[iy][ix].shbl1np.qu;
  bodygrid[iy][ix].shbl1np.qv=pls*bodygrid[iy][ix].shbl1np.qv;

  bodygrid[iy][ix].shbl1np.te=(bodygrid[iy][ix].shbl1np.te
  -sgm*bodygrid[iy][ix].jbnp.p)/(1.0+sgm);

  shu=bodygrid[iy][ix].shbl1np;

  shblToJbbl(&shu, &jbu);
  }
  }
  }
  }
  -------------------------------------------------------------------*/


