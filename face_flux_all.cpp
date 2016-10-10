
#include"non_uniform_grid.h"
#include"vec2d.h"
#include"vec4d.h"
extern OctCell *bodygrid;
extern double timesum;
extern double dt;
extern vector<Face> faces_comp;

//如果尖后缘在一个网格中,退化的网格那这断尖后缘被截断,网格当着计算单元
//多值点注意这个尖后缘点还有前缘点
//需要找到那个前缘点,尖后援所在的网格,那个网格特殊处理.有一个face的ghost cell需要平均上下面ghost cell的结果
//这个网格可能会标记为界面网格,或者就有可能就标记为流体网格,现在尖后缘被截断,那这个网格就是计算单元
//如果面两端有一个是需要处理的固体内部点,那所对应的ghost点与对称点的对称线以离按面扫描的面中点最近的那个翼形表面点的两端的两个折线断的圆作为径向对称线
//边界法线的处理这里借鉴dg的一篇文章关于高阶格式低阶边界近似的处理方法
//关于多值点最后处理通过适当的标记指比如标志为4来特殊处理

//x,y belong to [-1,1]X[-1,1] 需要将  [x_{i-1/2} x_{i+1/2}] => [-1,1]
double approx_sol(double x, double y, double *doftmp, int n=nDOF );
//直接x,y是原始计算区域的坐标 [x_{i-1/2} x_{i+1/2}] 
double n_approx_sol(double x,double y, double xc1, double yc1, double dx, double dy, double *doftmp, int n=nDOF );
double n_approx_sol(const double xi[2],  const double xc[2], double dx, double dy, double *doftmp, int n=nDOF );

void HLLC_flux(JBBL *ujbL, JBBL *ujbR, PXYZ *nml, double flxface[4]);
void HLLC_flux(SHBL &wl, SHBL &wr, PXYZ &normal, double flxface[4]);
void GLF_flux(SHBL &wl, SHBL &wr, PXYZ &nml, double flxface[4]); 
void getGhostValueMagazine(OctCell * ghostcell,int);
void getGhostValueRiemann(OctCell * ghostcell);

void face_flux(Face &fc)
{
     double uxl[4][3]={0.0},uxr[4][3]={0.0};
     
     OctCell *lsbl, *lsbl0;
     double xc0[2],xc1[2];
     
     double xii0[3][2]={0.0};
          
     double dx0,dy0,dx1,dy1;
     double dofp[4][6]={0.0};
     double dofn[4][6]={0.0};
     bool pflag2,nflag2; 
     lsbl0=fc.parent;
     lsbl=fc.neighbor;
     xc0[0]=lsbl0->xc1;
     xc0[1]=lsbl0->yc1;
     
     xc1[0]=lsbl->xc1;
     xc1[1]=lsbl->yc1;
     
     dx0=lsbl0->dx;
     dy0=lsbl0->dy;
     dx1=lsbl->dx;
     dy1=lsbl->dy;
      
     for(short j=0; j<3;++j){
        xii0[j][0]=fc.gs_pt[j][0];
        xii0[j][1]=fc.gs_pt[j][1];                 
      }

     if((!lsbl0->mvflg)||(!lsbl->mvflg)){
         for(short i1=0;i1<4;++i1){
             for(short j1=0;j1<6;++j1){
                 dofp[i1][j1]=lsbl0->dof[i1][j1];
                 dofn[i1][j1]=lsbl->dof[i1][j1];
             }
         }
     }
     else{
         pflag2=(lsbl0->flag==2);
         nflag2=(lsbl->flag==2);
         if((pflag2&& !nflag2) || (!pflag2&& nflag2)){
             //.....
             if(b_ARS){
                 if(pflag2) getGhostValueRiemann(lsbl0);
                 if(nflag2) getGhostValueRiemann(lsbl);
             }
             else{
                 if(pflag2) getGhostValueMagazine(lsbl0,mv_qvlvModify);
                 if(nflag2) getGhostValueMagazine(lsbl,mv_qvlvModify);
             }
             for(short i1=0;i1<4;++i1){
                 for(short j1=0;j1<6;++j1){
                     dofp[i1][j1]=lsbl0->dof[i1][j1];
                     dofn[i1][j1]=lsbl->dof[i1][j1];
                 }
             }
         }
         else {
             for(short i1=0;i1<4;++i1){
                 for(short j1=0;j1<6;++j1){
                     dofp[i1][j1]=lsbl0->dof[i1][j1];
                     dofn[i1][j1]=lsbl->dof[i1][j1];
                 }
             }
         }
     }

     /*
        for(short i1=0;i1<4;++i1){
        for(short j1=0;j1<6;++j1){
        dofp[i1][j1]=lsbl0->dof[i1][j1];
        dofn[i1][j1]=lsbl->dof[i1][j1];
        }
        }
        */
     for(short i=0; i< 4; i++)  {
         for(short j=0;j<3;j++) {
	    uxl[i][j]= n_approx_sol(xii0[j], xc0, dx0, dy0, dofp[i], nDOF );
	    uxr[i][j]= n_approx_sol(xii0[j], xc1, dx1, dy1, dofn[i], nDOF );
	    }
     }     
     
     SHBL ushl[3], ushr[3];
   //  JBBL ujbl[3], ujbr[3];     
     
     for(short j=0; j<3;j++) {
          ushl[j].q=uxl[0][j];
          ushl[j].qu=uxl[1][j];
          ushl[j].qv=uxl[2][j];         
          ushl[j].te=uxl[3][j];
          
          ushr[j].q=uxr[0][j];
          ushr[j].qu=uxr[1][j];
          ushr[j].qv=uxr[2][j];         
          ushr[j].te=uxr[3][j];
      }     
     
     PXYZ normal;
     
     normal.x=fc.nml[0];
     normal.y=fc.nml[1];
     
     double flx[3][4]={0.0};
     for(short j=0; j<3;++j) {
      //    shblToJbbl(&ushl[j], &ujbl[j]);
      //    shblToJbbl(&ushr[j], &ujbr[j]); 
  //        HLLC_flux(ujbl+j, ujbr+j, &normal, flx[j]); 

  //void HLLC_flux(SHBL &wl, SHBL &wr, PXYZ &normal, double flxface[4])  
         if(HLLCFLUX==0)
             GLF_flux(ushl[j], ushr[j], normal, flx[j]); 
         else
             HLLC_flux(ushl[j], ushr[j], normal, flx[j]); 
     }
     
     for(short j=0;j<3; ++j) {
        for(short i=0; i<4; ++i) {
           fc.faceflux[j][i]=flx[j][i];
         }
     }
/*
#ifdef KXRCF_LIM
     double u1,v1,u0,v0;
     double rtmp;
     JBBL uljb[3],urjb[3];
     double prel[3], prer[3];

     for(short j=0; j<3;++j) {
          shblToJbbl(&ushl[j], &uljb[j]);
          shblToJbbl(&ushr[j], &urjb[j]); 
          prel[j]=uljb[j].p;
          prer[j]=urjb[j].p;
     }
     u0=dofp[1][0]/dofp[0][0];
     v0=dofp[2][0]/dofp[0][0];
     u1=dofn[1][0]/dofn[0][0];
     v1=dofn[2][0]/dofn[0][0];
//den
 //    rtmp=(wei[0]*(uxl[0][0]-uxr[0][0])+wei[1]*(uxl[0][1]-uxr[0][1])
   //      +wei[2]*(uxl[0][2]-uxr[0][2]))*fc.area*0.5;
     //    
     //entropy
     rtmp=(wei[0]*(prel[0]/pow(uxl[0][0],GAMMA)-prer[0]/pow(uxr[0][0],GAMMA))
         +wei[1]*(prel[1]/pow(uxl[0][1],GAMMA)-prer[1]/pow(uxr[0][1],GAMMA)) 
         +wei[2]*(prel[2]/pow(uxl[0][2],GAMMA)-prer[2]/pow(uxr[0][2],GAMMA)) )*fc.area*0.5;
     if(u0*normal.x+v0*normal.y<0.0){
         lsbl0->kxq+=rtmp;
         lsbl0->inlen+=fc.area;
     }

     if(u1*normal.x+v1*normal.y>0.0){
         lsbl->kxq-=rtmp;
         lsbl->inlen+=fc.area;
     }

#endif
*/

 }      
  
void face_kxrcf(Face &fc)
{
     double uxl[4][3]={0.0},uxr[4][3]={0.0};
     
     OctCell *lsbl, *lsbl0;
     double xc0[2],xc1[2];
     
     double xii0[3][2]={0.0};
          
     double dx0,dy0,dx1,dy1;
     double dofp[4][6]={0.0};
     double dofn[4][6]={0.0};
     bool pflag2,nflag2; 
     lsbl0=fc.parent;
     lsbl=fc.neighbor;
     xc0[0]=lsbl0->xc1;
     xc0[1]=lsbl0->yc1;
     
     xc1[0]=lsbl->xc1;
     xc1[1]=lsbl->yc1;
     
     dx0=lsbl0->dx;
     dy0=lsbl0->dy;
     dx1=lsbl->dx;
     dy1=lsbl->dy;
      
     for(short j=0; j<3;++j){
        xii0[j][0]=fc.gs_pt[j][0];
        xii0[j][1]=fc.gs_pt[j][1];                 
      }

     if((!lsbl0->mvflg)||(!lsbl->mvflg)){
         for(short i1=0;i1<4;++i1){
             for(short j1=0;j1<6;++j1){
                 dofp[i1][j1]=lsbl0->dof[i1][j1];
                 dofn[i1][j1]=lsbl->dof[i1][j1];
             }
         }
     }
     else{
         pflag2=(lsbl0->flag==2);
         nflag2=(lsbl->flag==2);
         if((pflag2&& !nflag2) || (!pflag2&& nflag2)){
             //.....
             if(pflag2) getGhostValueMagazine(lsbl0,mv_qvlvModify);
             if(nflag2) getGhostValueMagazine(lsbl,mv_qvlvModify);
             for(short i1=0;i1<4;++i1){
                 for(short j1=0;j1<6;++j1){
                     dofp[i1][j1]=lsbl0->dof[i1][j1];
                     dofn[i1][j1]=lsbl->dof[i1][j1];
                 }
             }
         }
         else {
             for(short i1=0;i1<4;++i1){
                 for(short j1=0;j1<6;++j1){
                     dofp[i1][j1]=lsbl0->dof[i1][j1];
                     dofn[i1][j1]=lsbl->dof[i1][j1];
                 }
             }
         }
     }

     /*
        for(short i1=0;i1<4;++i1){
        for(short j1=0;j1<6;++j1){
        dofp[i1][j1]=lsbl0->dof[i1][j1];
        dofn[i1][j1]=lsbl->dof[i1][j1];
        }
        }
        */
     for(short i=0; i< 4; i++)  {
         for(short j=0;j<3;j++) {
	    uxl[i][j]= n_approx_sol(xii0[j], xc0, dx0, dy0, dofp[i], nDOF );
	    uxr[i][j]= n_approx_sol(xii0[j], xc1, dx1, dy1, dofn[i], nDOF );
	    }
     }     
     
     SHBL ushl[3], ushr[3];
   //  JBBL ujbl[3], ujbr[3];     
     
     for(short j=0; j<3;j++) {
          ushl[j].q=uxl[0][j];
          ushl[j].qu=uxl[1][j];
          ushl[j].qv=uxl[2][j];         
          ushl[j].te=uxl[3][j];
          
          ushr[j].q=uxr[0][j];
          ushr[j].qu=uxr[1][j];
          ushr[j].qv=uxr[2][j];         
          ushr[j].te=uxr[3][j];
      }     
     
     PXYZ normal;
     
     normal.x=fc.nml[0];
     normal.y=fc.nml[1];
     
     double u1,v1,u0,v0;
     double rtmp;
     JBBL uljb[3],urjb[3];
     double prel[3], prer[3];

     for(short j=0; j<3;++j) {
          shblToJbbl(&ushl[j], &uljb[j]);
          shblToJbbl(&ushr[j], &urjb[j]); 
          prel[j]=uljb[j].p;
          prer[j]=urjb[j].p;
     }
     u0=dofp[1][0]/dofp[0][0];
     v0=dofp[2][0]/dofp[0][0];
     u1=dofn[1][0]/dofn[0][0];
     v1=dofn[2][0]/dofn[0][0];
/*den
     rtmp=(wei[0]*(uxl[0][0]-uxr[0][0])+wei[1]*(uxl[0][1]-uxr[0][1])
         +wei[2]*(uxl[0][2]-uxr[0][2]))*fc.area*0.5;
         */
     //entropy
     rtmp=(wei[0]*(prel[0]/pow(uxl[0][0],GAMMA)-prer[0]/pow(uxr[0][0],GAMMA))
         +wei[1]*(prel[1]/pow(uxl[0][1],GAMMA)-prer[1]/pow(uxr[0][1],GAMMA)) 
         +wei[2]*(prel[2]/pow(uxl[0][2],GAMMA)-prer[2]/pow(uxr[0][2],GAMMA)) )*fc.area*0.5;
     if(u0*normal.x+v0*normal.y<0.0){
         lsbl0->kxq+=rtmp;
         lsbl0->inlen+=fc.area;
     }

     if(u1*normal.x+v1*normal.y>0.0){
         lsbl->kxq-=rtmp;
         lsbl->inlen+=fc.area;
     }

 }      
 
void face_kxrcf_all( )
{
   for(vector<Face>::size_type isz=0; isz!=faces_comp.size(); ++isz){
       face_kxrcf(faces_comp[isz]);
   }   
}

void set_max_cv();
void face_flux_all( )
{
    if(HLLCFLUX==0)
        set_max_cv();
    for(vector<Face>::size_type isz=0; isz!=faces_comp.size(); ++isz){
        face_flux(faces_comp[isz]);
    }   
}


extern Node *HeadListAllGrid;
void kxrcf_relate()
{
    Node *current;
    current = HeadListAllGrid;  
    OctCell *pcell0;
    double rtp, rtm;
    double uapp[4][3][3];
    double pre[3][3]={0.0};
    while(current != NULL)
    {		
        pcell0 = current->cell;
        if(pcell0->flag  == 0 && current->flg<=2 ) {
            pcell0->kxq=0.0;
            pcell0->inlen=0.0;
            rtp=1.0e-13;
            for(short j=0; j<4; j++) {
                for(short imm=0; imm<3;imm++) {
                    for(short jmm=0; jmm<3; jmm++) {
                        uapp[j][jmm][imm]=approx_sol(gspt[jmm],
                            gspt[imm], pcell0->dof[j], nDOF);  //u(x,y)
                       // if(rtp<uapp[jmm][imm]) rtp=uapp[0][jmm][imm];
                       // density
                    }
                }
            }
            //entropy
            for(short imm=0; imm<3;imm++) {
                for(short jmm=0; jmm<3; jmm++) {
                    pre[jmm][imm]=GAM11 *(uapp[3][jmm][imm]-0.5*(uapp[1][jmm][imm]*uapp[1][jmm][imm]/uapp[0][jmm][imm]
                            + uapp[2][jmm][imm]*uapp[2][jmm][imm]/uapp[0][jmm][imm]));	          
                    rtm=pre[jmm][imm]/pow(uapp[0][jmm][imm],GAMMA);
                    if(rtp<rtm) rtp=rtm;
                }
            }
            pcell0->maxqp=rtp;
        }
        current = current->next;
    }
}

extern double value_kxrcf_amr;
void set_kxrcf()
{
    Node *current;
    current = HeadListAllGrid;  
    OctCell *pcell0;
    double rtp;
  //  int pol=2;
  //  pol=(nDOF==6?2:1);
    while(current != NULL)
    {		
        pcell0 = current->cell;
        if(pcell0->flag  == 0 && current->flg<=2 ) {
          //  rtp=max2(pcell0->dx,pcell0->dy);
            rtp=0.5*sqrt(pcell0->dx*pcell0->dx
                +pcell0->dy*pcell0->dy);
            pcell0->kxrcf=fabs(pcell0->kxq)/pow(rtp,1.5)/
            pcell0->inlen/pcell0->maxqp;
            if(pcell0->kxrcf>1.0) pcell0->trb=true;
            else pcell0->trb=false;
            if(pcell0->kxrcf>value_kxrcf_amr) pcell0->trbamr=true;
            else pcell0->trbamr=false;
        }
        current = current->next;
    }
}



