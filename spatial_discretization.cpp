

#include"non_uniform_grid.h"
#include"findneighbor.h"
extern OctCell *bodygrid;
extern Node *HeadListAllGrid;
//x,y belong to [-1,1]X[-1,1] 需要将  [x_{i-1/2} x_{i+1/2}] => [-1,1]
double approx_sol(double x, double y, double *doftmp, int n=nDOF );
//直接x,y是原始计算区域的坐标 [x_{i-1/2} x_{i+1/2}] 
double n_approx_sol(double x,double y, double xc1, double yc1, double dx, double dy, double *doftmp, int n=nDOF );
double n_approx_sol(const double xi[2],  const double xc[2], double dx, double dy, double *doftmp, int n=nDOF );

extern vector<Face> faces_comp; 

void spatial_discretize_f(Face &fc)
{
   OctCell *un;
   OctCell *un1;
 //  OctCell *unp1x, *unp1y, *unm1x, *unm1y;

//   double tmp92=0.0;
//   double tmp93=0.0;
//   double tmp94=0.0;  
//   double tmpx=0.0, tmpy=0.0;
   
   un=fc.parent;
   un1=fc.neighbor;
   double xc1,yc1, xc2,yc2;
   double dx1,dy1,dy2,dx2;     
//   double uapp[4][3][3]={0.0};
//   double uapp1[4][3][3]={0.0};
    
   xc1=un->xc1;
   yc1=un->yc1;
   dx1=un->dx;
   dy1=un->dy;
   
   xc2=un1->xc1;
   yc2=un1->yc1;
   dx2=un1->dx;
   dy2=un1->dy;  
         
  double fflx[4][6]={0.0};
  double rfflx[4][6]={0.0};
			      	  
   for (short j=0; j<4; j++) {
       fflx[j][0]=(wei[0]*fc.faceflux[0][j]+wei[1]*fc.faceflux[1][j]+wei[2]*fc.faceflux[2][j])*fc.area*0.5;
       fflx[j][1]=(wei[0]*fc.faceflux[0][j]*phix(fc.gs_pt[0][0],xc1,dx1)+wei[1]*fc.faceflux[1][j]*phix(fc.gs_pt[1][0],xc1,dx1)
             +wei[2]*fc.faceflux[2][j]*phix(fc.gs_pt[2][0],xc1,dx1))*fc.area*0.5;       
        
       fflx[j][2]=(wei[0]*fc.faceflux[0][j]*psiy(fc.gs_pt[0][1],yc1,dy1)+wei[1]*fc.faceflux[1][j]*psiy(fc.gs_pt[1][1],yc1,dy1)
             +wei[2]*fc.faceflux[2][j]*psiy(fc.gs_pt[2][1],yc1,dy1))*fc.area*0.5;         
        
//    if(nDOF==6) {        
       fflx[j][3]=(wei[0]*fc.faceflux[0][j]*phix_psiy(fc.gs_pt[0][0], fc.gs_pt[0][1],xc1,yc1,dx1, dy1)
             +wei[1]*fc.faceflux[1][j]*phix_psiy(fc.gs_pt[1][0], fc.gs_pt[1][1],xc1,yc1,dx1,dy1)
             +wei[2]*fc.faceflux[2][j]*phix_psiy(fc.gs_pt[2][0],fc.gs_pt[2][1],xc1,yc1,dx1, dy1))*fc.area*0.5;         
             
       fflx[j][4]=(wei[0]*fc.faceflux[0][j]*phix2m(fc.gs_pt[0][0],xc1,dx1)+wei[1]*fc.faceflux[1][j]*phix2m(fc.gs_pt[1][0],xc1,dx1)
             +wei[2]*fc.faceflux[2][j]*phix2m(fc.gs_pt[2][0],xc1,dx1))*fc.area*0.5;       
          
       fflx[j][5]=(wei[0]*fc.faceflux[0][j]*psiy2m(fc.gs_pt[0][1],yc1,dy1)+wei[1]*fc.faceflux[1][j]*psiy2m(fc.gs_pt[1][1],yc1,dy1)
             +wei[2]*fc.faceflux[2][j]*psiy2m(fc.gs_pt[2][1],yc1,dy1))*fc.area*0.5;   
  //   }
//right face   
       rfflx[j][0]=fflx[j][0];
       rfflx[j][1]=(wei[0]*fc.faceflux[0][j]*phix(fc.gs_pt[0][0],xc2,dx2)+wei[1]*fc.faceflux[1][j]*phix(fc.gs_pt[1][0],xc2,dx2)
             +wei[2]*fc.faceflux[2][j]*phix(fc.gs_pt[2][0],xc2,dx2))*fc.area*0.5;       
        
       rfflx[j][2]=(wei[0]*fc.faceflux[0][j]*psiy(fc.gs_pt[0][1],yc2,dy2)+wei[1]*fc.faceflux[1][j]*psiy(fc.gs_pt[1][1],yc2,dy2)
             +wei[2]*fc.faceflux[2][j]*psiy(fc.gs_pt[2][1],yc2,dy2))*fc.area*0.5;         
        
//    if(nDOF==6) {        
       rfflx[j][3]=(wei[0]*fc.faceflux[0][j]*phix_psiy(fc.gs_pt[0][0], fc.gs_pt[0][1],xc2,yc2,dx2, dy2)
             +wei[1]*fc.faceflux[1][j]*phix_psiy(fc.gs_pt[1][0], fc.gs_pt[1][1],xc2,yc2,dx2,dy2)
             +wei[2]*fc.faceflux[2][j]*phix_psiy(fc.gs_pt[2][0],fc.gs_pt[2][1],xc2,yc2,dx2, dy2))*fc.area*0.5;         
             
       rfflx[j][4]=(wei[0]*fc.faceflux[0][j]*phix2m(fc.gs_pt[0][0],xc2,dx2)+wei[1]*fc.faceflux[1][j]*phix2m(fc.gs_pt[1][0],xc2,dx2)
             +wei[2]*fc.faceflux[2][j]*phix2m(fc.gs_pt[2][0],xc2,dx2))*fc.area*0.5;       
          
       rfflx[j][5]=(wei[0]*fc.faceflux[0][j]*psiy2m(fc.gs_pt[0][1],yc2,dy2)+wei[1]*fc.faceflux[1][j]*psiy2m(fc.gs_pt[1][1],yc2,dy2)
             +wei[2]*fc.faceflux[2][j]*psiy2m(fc.gs_pt[2][1],yc2,dy2))*fc.area*0.5;   
  //   }  
  }          
   	             
   //face flux
   for(short j=0;j<4;++j){
      for(short imm=0; imm<6;++imm) {
         un->residual[j][imm] += fflx[j][imm];   
         un1->residual[j][imm] -= rfflx[j][imm]; 
      }
  }	    

//end face flux
}


void spatial_discretize_v(OctCell* un)
{
   double tmp92=0.0;
   double tmp93=0.0;
   double tmp94=0.0;  
   double tmpx=0.0, tmpy=0.0;
   
   double uapp[4][3][3]={0.0};
   double dx,dy, dxpdy;
   dx=un->dx;
   dy=un->dy;
   dxpdy=dx*dy;
   
   for(short j=0; j<4; j++) {
      for(short imm=0; imm<3;imm++) {
	   for(short jmm=0; jmm<3; jmm++) {
	       uapp[j][jmm][imm]=approx_sol(gspt[jmm], gspt[imm], un->dof[j], nDOF);  //u(x,y)
	    }
	}
   }
	
  double pre[3][3]={0.0};
  double pre1[3][3]={0.0};	
  	
  for(short imm=0; imm<3;imm++) {
     for(short jmm=0; jmm<3; jmm++) {
	      pre[jmm][imm]=GAM11 *(uapp[3][jmm][imm]-0.5*(uapp[1][jmm][imm]*uapp[1][jmm][imm]/uapp[0][jmm][imm]
	          + uapp[2][jmm][imm]*uapp[2][jmm][imm]/uapp[0][jmm][imm]));	          
     }
  }
          
  double flxfunx[4][3][3]={0.0};
  double flxfuny[4][3][3]={0.0};
  
  for(short imm=0; imm<3;imm++) {
     for(short jmm=0; jmm<3; jmm++) {

	         flxfunx[0][jmm][imm]=uapp[1][jmm][imm]; 
	         flxfuny[0][jmm][imm]=uapp[2][jmm][imm]; 	
	                  
	         flxfunx[1][jmm][imm]=uapp[1][jmm][imm]*uapp[1][jmm][imm]/uapp[0][jmm][imm] + pre[jmm][imm];

	         flxfunx[2][jmm][imm] = flxfuny[1][jmm][imm] = uapp[1][jmm][imm]*uapp[2][jmm][imm]/uapp[0][jmm][imm];
	         
	         flxfuny[2][jmm][imm]=uapp[2][jmm][imm]*uapp[2][jmm][imm]/uapp[0][jmm][imm] + pre[jmm][imm];	
	         
	         flxfunx[3][jmm][imm] = (uapp[3][jmm][imm] + pre[jmm][imm])* uapp[1][jmm][imm]/uapp[0][jmm][imm];
	         flxfuny[3][jmm][imm] = (uapp[3][jmm][imm] + pre[jmm][imm])* uapp[2][jmm][imm]/uapp[0][jmm][imm];
	           	         	              	     
	      }
  }
  
//体积积分项
     
// euler equation for aerodynamics 
//uapp[j][jmm][imm]
    for(short j=0; j<4;++j) {
       tmp92=0.0;
	tmp93=0.0;
	tmp94=0.0;	
	tmpx=tmpy=0.0;	
	
//base function: phi_i(x)  &  psi_j(y)
	for(short imm=0; imm < 3;imm++) {
	   for(short jmm=0; jmm < 3; jmm++) {
	      tmp92 += flxfunx[j][jmm][imm]*wei[imm]*wei[jmm];
	      
	      tmp93 += flxfuny[j][jmm][imm]*wei[imm]*wei[jmm];
	      }
	 }	 

	 tmp92 *=  0.25*dxpdy;
	 tmp93 *=  0.25*dxpdy;	 
	       
        un->residual[j][1] -= 2.0/dx * tmp92;
        un->residual[j][2] -= 2.0/dy * tmp93;
        
        tmp92=tmp93=0.0;
  
  //base function: phi_i(x)*psi_j(y), phi_i(x)^2-1.0/3.0, &  psi_j(y)^2-1.0/3
     if(nDOF==6) {
	 for(short imm=0; imm<3;imm++) {
	   for(short jmm=0; jmm<3; jmm++) {
	      tmp93 += flxfuny[j][jmm][imm] *gspt[jmm]*wei[imm]*wei[jmm];   //fy*x 
	      
	      tmp94 += flxfunx[j][jmm][imm] *gspt[imm]*wei[imm]*wei[jmm];   //fx*y 
	      
	      tmpx += flxfunx[j][jmm][imm] *gspt[jmm]*wei[imm]*wei[jmm];  //fx*x
	      tmpy += flxfuny[j][jmm][imm] *gspt[imm]*wei[imm]*wei[jmm];	//fy*y      
	      }
	  }

	  tmp94 *= 0.25*dxpdy;
	  tmp93 *= 0.25*dxpdy;
	  tmpx *= 0.25*dxpdy;
	  tmpy *= 0.25*dxpdy;

	  un->residual[j][3] -= 2.0/dx * tmp94 + 2.0/dy * tmp93;
	
	  un->residual[j][4] -= 4.0/dx * tmpx;
	  un->residual[j][5] -= 4.0/dy * tmpy;  
	  	
       }
    } //end for
//=====
    
}          

void spatial_discretization_all()
{
    Node *current;
    current = HeadListAllGrid;  
    OctCell *pcell0;
    while(current != NULL)
    {		
        pcell0 = current->cell;
        if(pcell0->flag % 2 == 0 && current->flg<=2 ) {
            for(int im=0;im<4; ++im){
                for(int jm=0;jm<6; ++jm){
                    pcell0->residual[im][jm]=0.0;
                }
            }
        }

        current = current->next;
    }


    for(vector<Face>::size_type isz=0; isz!=faces_comp.size(); ++isz){
        spatial_discretize_f(faces_comp[isz]);
    }

    current = HeadListAllGrid;     
    while(current != NULL)
    {		
        pcell0 = current->cell;
        if(pcell0->flag == 0 && current->flg<=2 ) 
            spatial_discretize_v(pcell0);
        current = current->next;
    }
}

#ifdef RES_SMOOTH
const double sm_ep=500.1;  //0.5-0.8
void resd_smoothing()
{
    OctCell *pc0,*pc1;
    Node *current;
    current = HeadListAllGrid;  
    OctCell *pcell0;
    double vol0,vol1;
    while(current != NULL)
    {		
        pcell0 = current->cell;
        if(pcell0->flag  == 0 && current->flg<=2 ) {
            for(int im=0;im<4; ++im){
                for(int jm=0;jm<6; ++jm){
                    pcell0->resd1[im][jm]=pcell0->residual[im][jm];

             //       pcell0->resd0[im][jm]=0.0;

                }
            }
        }

        current = current->next;
    }

    OctCell *pcs,*pcn;
    double renb[4][6],renb1[4][6];
    double renbs[4][6],renbn[4][6];
    double rtmx,rtmy;
    for(int itm=0;itm<2;++itm){
        current = HeadListAllGrid;  
        while(current != NULL)
        {		
            pcell0 = current->cell;
            if(pcell0->flag  == 0 && current->flg<=2 ) {
                pc0=WestNeighbor(pcell0);
                pc1=EastNeighbor(pcell0);
                pcs=SouthNeighbor(pcell0);
                pcn=NorthNeighbor(pcell0);

                rtmx=rtmy=0.0;
                if(pc0->reflag==1){
                    for(int iu=0;iu<4;++iu)
                        for(int ju=0;ju<6;++ju)
                            renb[iu][ju]=0.0;
                }
                else {
                    ++rtmx;
                    for(int iu=0;iu<4;++iu)
                        for(int ju=0;ju<6;++ju)
                            renb[iu][ju]=pc0->resd1[iu][ju];
                }

                if(pc1->reflag==1){
                    for(int iu=0;iu<4;++iu)
                        for(int ju=0;ju<6;++ju)
                            renb1[iu][ju]=0.0;
                }
                else {
                    ++rtmx;
                    for(int iu=0;iu<4;++iu)
                        for(int ju=0;ju<6;++ju)
                            renb1[iu][ju]=pc1->resd1[iu][ju];
                }


                if(pcs->reflag==1){
                    for(int iu=0;iu<4;++iu)
                        for(int ju=0;ju<6;++ju)
                            renbs[iu][ju]=0.0;
                }
                else {
                    ++rtmy;
                    for(int iu=0;iu<4;++iu)
                        for(int ju=0;ju<6;++ju)
                            renbs[iu][ju]=pcs->resd1[iu][ju];
                }

                if(pcn->reflag==1){
                    for(int iu=0;iu<4;++iu)
                        for(int ju=0;ju<6;++ju)
                            renbn[iu][ju]=0.0;
                }
                else {
                    ++rtmy;
                    for(int iu=0;iu<4;++iu)
                        for(int ju=0;ju<6;++ju)
                            renbn[iu][ju]=pcn->resd1[iu][ju];
                }

                for(int im=0;im<4; ++im){
                    for(int jm=0;jm<6; ++jm){
                        pcell0->resd1[im][jm]
                            =(pcell0->resd1[im][jm]+sm_ep*(renb[im][jm]+renb1[im][jm]))
                            /(1.0+sm_ep*rtmx);
                    }
                }

                for(int im=0;im<4; ++im){
                    for(int jm=0;jm<6; ++jm){
                        pcell0->resd1[im][jm]
                            =(pcell0->resd1[im][jm]+sm_ep*(renbs[im][jm]+renbn[im][jm]))
                            /(1.0+sm_ep*rtmy);
                    }
                }
            }

            current = current->next;
        }
    }

/*
    for(int itmp=0;itmp<2;++itmp){
        for(vector<Face>::size_type isz=0; isz!=faces_comp.size(); ++isz){
            pc0=faces_comp[isz].parent;
            vol0=pc0->dx*pc0->dy;
            vol0=1.;
            pc1=faces_comp[isz].neighbor;
            vol1=pc1->dx*pc1->dy;
            vol1=1.;
            for(int im=0;im<4;++im){
                for(int jm=0;jm<6;++jm){
                    if(pc0->flag==0){
                        pc0->resd0[im][jm]+=vol0/vol1*pc1->resd1[im][jm]*sm_ep;
                    }
                    if(pc1->flag==0){
                        pc1->resd0[im][jm]+=vol1/vol0*pc0->resd1[im][jm]*sm_ep;
                    }
                }
            }
        }

        current = HeadListAllGrid;  
        while(current != NULL)
        {		
            pcell0 = current->cell;
            if(pcell0->flag  == 0 && current->flg<=2 ) {
                for(int im=0;im<4; ++im){
                    for(int jm=0;jm<6; ++jm){
                        pcell0->resd1[im][jm]
                            =(pcell0->resd1[im][jm]+pcell0->resd0[im][jm])
                            /(1.0+sm_ep*pcell0->p_node->n_nb);
                    }
                }
            }

            current = current->next;
        }

    */
        current = HeadListAllGrid;  
        while(current != NULL)
        {		
            pcell0 = current->cell;
            if(pcell0->flag  == 0 && current->flg<=2 ) {
  //              for(int im=0;im<4; ++im){
//                    for(int jm=0;jm<6; ++jm){
                    //    if(fabs(pcell0->resd1[im][jm])
                      //      <fabs(pcell0->residual[im][jm]))
    //                        pcell0->residual[im][jm]=pcell0->resd1[im][jm];
      //              }
        //        }

                for(int im=0;im<4; ++im){
                    if(fabs(pcell0->resd1[im][0])
                        <fabs(pcell0->residual[im][0]))
                        pcell0->residual[im][0]=pcell0->resd1[im][0];
                }

            }

            current = current->next;
        }
    //}
}
#endif


