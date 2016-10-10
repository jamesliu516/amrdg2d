#include"non_uniform_grid.h"

extern OctCell *bodygrid;

//extern REAL initQ, uFree, vFree, initP,initK,initOmg;
extern Node *HeadListAllGrid;
extern JBBL State_l;
extern JBBL State_r;
void jbblToShbl(JBBL *, SHBL *);
void left_shock_sol(double qr,double ur, double pr,double Ms,
	double &ql, double &ul, double &pl)//Ms shock wave Mach number
{
	double Mr;
	double ar;
	double gr;
	double mrs2;
	gr=GAMMA;
	ar=sqrt(gr*pr/qr);
	Mr=ur/ar;
	mrs2=(Mr-Ms)*(Mr-Ms);
	ql=qr*(gr+1)*mrs2/((gr-1)*mrs2+2.0);
	pl=pr*(2.0*gr*mrs2-(gr-1))/(gr+1);
	ul=(1.0-qr/ql)*Ms*ar+ur*qr/ql;
}
void right_shock_sol(double ql,double ul, double pl,double Ms,
	double &qr, double &ur, double &pr)//Ms shock wave Mach number
{
	double Ml;
	double al;
	double gr;
	double mrs2;
	gr=GAMMA;
	al=sqrt(gr*pl/ql);
	Ml=ul/al;
	mrs2=(Ml-Ms)*(Ml-Ms);
	qr=ql*(gr+1)*mrs2/((gr-1)*mrs2+2.0);
	pr=pl*(2.0*gr*mrs2-(gr-1))/(gr+1);
	ur=(1.0-ql/qr)*Ms*al+ul*ql/qr;
}
	
//(x,y) belong to [-1,1]x[-1,1] 在正规标准正方形元上定义
double u0xy(double x, double y, double xc, double yc, double dx, double dy, int stp)
{
    double r;

    double x0, y0, hxt0;
    x0=0.5*x*dx+xc;
    y0=0.5*y*dy+yc; 

    JBBL ujb[4];
    SHBL ush[4];
   
    JBBL jbu1;
    SHBL shu1;
/*
//    unsteady shock vortex interaction 
	double ur,pr,qr,ul,pl,ql;
    State_l.q=1.0;
    State_l.p=initP;
    State_l.u=0.0;
    State_l.v=0.0;

	qr=jbu1.q=State_r.q;
	pr=jbu1.p=State_r.p;
	ur=jbu1.u=jbu1.v=0.0;

	jbblToShbl(&jbu1, &shu1);
    JBBL jbu0;
    SHBL shu0;
	left_shock_sol(qr,ur,pr,FreeMa,ql,ul,pl);
    State_l.q=ql;
    State_l.p=pl;
	jbu0.q=ql;
	State_l.u=jbu0.u=ul;
	State_l.v=jbu0.v=0.0;
	jbu0.p=pl;
	jbblToShbl(&jbu0, &shu0);
	switch(stp){
	case 1:
		if(x0<0.385) r=shu0.q;
		else r=shu1.q;
		break;
	case 2:
		if(x0<0.385)  r=shu0.qu;
		else  r=shu1.qu;
		break;
	case 3:
		if(x0<0.385)  r=shu0.qv;
		else  r=shu1.qv;
		break;
	case 4:
		if(x0<0.385)  r=shu0.te;
		else  r=shu1.te;
		break;
	}
    */
    //-----------------------------------
/*
//    unsteady shock cylinder
	double ur,pr,qr,ul,pl,ql;
    State_r.q=1.0;
    State_r.p=initP;
    State_r.u=0.0;
    State_r.v=0.0;

	qr=jbu1.q=State_r.q;
	pr=jbu1.p=State_r.p;
	ur=jbu1.u=jbu1.v=0.0;

	jbblToShbl(&jbu1, &shu1);
    JBBL jbu0;
    SHBL shu0;
	left_shock_sol(qr,ur,pr,FreeMa,ql,ul,pl);
    State_l.q=ql;
    State_l.p=pl;
	jbu0.q=ql;
	State_l.u=jbu0.u=ul;
	State_l.v=jbu0.v=0.0;
	jbu0.p=pl;
	jbblToShbl(&jbu0, &shu0);
	switch(stp){
	case 1:
		if(x0<0.385) r=shu0.q;
		else r=shu1.q;
		break;
	case 2:
		if(x0<0.385)  r=shu0.qu;
		else  r=shu1.qu;
		break;
	case 3:
		if(x0<0.385)  r=shu0.qv;
		else  r=shu1.qv;
		break;
	case 4:
		if(x0<0.385)  r=shu0.te;
		else  r=shu1.te;
		break;
	}
    */

      //steady flow 
    jbu1.q=initQ; 
    jbu1.u=uFree;
    jbu1.v=vFree;

    jbu1.p=initP;
    jbblToShbl(&jbu1, &shu1);
        
   if(stp==1) {  

      r=shu1.q;
    }  
   else if(stp==2) {
      r=shu1.qu;

   }
   else if( stp==3) {
      r=shu1.qv;
   }
   else if(stp==4) {
      r=shu1.te;

   }     
    
   /*
// test the order
    JBBL jbu1;
    SHBL shu1;
 
    jbu1.q=1.0+0.2*sin(PI*(0.5*x*dx+xc+0.5*y*dy+yc)); 
    jbu1.u=0.7;
    jbu1.v=0.3;
    jbu1.p=1.0;
    jbblToShbl(&jbu1, &shu1);
        
   if(stp==1) {  
      r=shu1.q;
    }  
   else if(stp==2) {
      r=shu1.qu;
   }
   else if( stp==3) {
      r=shu1.qv;
   }

   else if(stp==4) {
      r=shu1.te;

   }     
    
      
    //2d Riemann
 
        ___________________
           |        |        |
           |        |        |
           |    1   |   3    |
           |________|________|
           |        |        |
           |    0   |    2   |
           |        |        |
           |________|________|
           

 
    ujb[0].q=1.1;
    ujb[0].u=0.8939;
    ujb[0].v=0.8939;
    ujb[0].p=1.1;
 
    ujb[1].q=0.5065;
    ujb[1].u=0.8939;
    ujb[1].v=0.0;
    ujb[1].p=0.350;
    
     ujb[2].q=0.5065;
    ujb[2].u=0.0;
    ujb[2].v=0.8939;
    ujb[2].p=0.350;    

    ujb[3].q=1.1;
    ujb[3].u=0.0;
    ujb[3].v=0.0;
    ujb[3].p=1.1; 
    */    
  /*   
    ujb[0].q=1.0;
    ujb[0].u=0.75;
    ujb[0].v=0.5;
    ujb[0].p=1.0;
 
    ujb[1].q=2.0;
    ujb[1].u=-0.75;
    ujb[1].v=0.5;
    ujb[1].p=1.0;
    
     ujb[2].q=3.0;
    ujb[2].u=0.75;
    ujb[2].v=-0.5;
    ujb[2].p=1.0;    

    ujb[3].q=1.0;
    ujb[3].u=-0.75;
    ujb[3].v=-0.5;
    ujb[3].p=1.0;     
 
   

    
   for( int itp=0; itp <4; itp++)
     jbblToShbl(ujb+itp, ush+itp);
    
   if(stp==1) {  
       if( x0<= 0.5 && y0 <= 0.5 ) { r=ush[0].q; } 
       else if(x0 <=0.5 && y0 >0.5)
       { r=ush[1].q;}
       else if(x0>0.5 && y0 <= 0.5)
       { r=ush[2].q;}
       else 
       { r=ush[3].q;}
    }  
   else if(stp==2) {
       if( x0<= 0.5 && y0 <= 0.5 ) { r=ush[0].qu; } 
       else if(x0 <=0.5 && y0 >0.5)
       { r=ush[1].qu;}
       else if(x0>0.5 && y0 <= 0.5)
       { r=ush[2].qu;}
       else 
       { r=ush[3].qu;}
   }
   else if( stp==3) {
       if( x0<= 0.5 && y0 <= 0.5 ) { r=ush[0].qv; } 
       else if(x0 <=0.5 && y0 >0.5)
       { r=ush[1].qv;}
       else if(x0>0.5 && y0 <= 0.5)
       { r=ush[2].qv;}
       else 
       { r=ush[3].qv;}
   }
   else if(stp==4) {
       if( x0<= 0.5 && y0 <= 0.5 ) { r=ush[0].te; } 
       else if(x0 <=0.5 && y0 >0.5)
       { r=ush[1].te;}
       else if(x0>0.5 && y0 <= 0.5)
       { r=ush[2].te;}
       else 
       { r=ush[3].te;}
   }  
    */      
   
 //double Mach reflection problem
   /*    
     hxt0=sqrt(3.0)*(x0-1.0/6.0);
    
     if(stp==1) {  
       if( y0-hxt0>= 0.0 || (y0<0.0 && x0<1.0/6.0) ) { r=8.0; } 
       else{ r=1.4;}
    }  
   else if(stp==2) {
      if(  y0-hxt0>= 0.0 || (y0<0.0 && x0<1.0/6.0))   r=57.1597;
      else
         r=0.0;
   }
   else if( stp==3) {
       if(  y0-hxt0>= 0.0 || (y0<0.0 && x0<1.0/6.0) )   r=-33.0012;
      else
         r=0.0;  
   }
   else if(stp==4) {
       if(  y0-hxt0>= 0.0 || (y0<0.0 && x0<1.0/6.0))   r=563.544;
      else
         r=2.5;
   } 

 
//Explosion Test in Two Space Dimensions see Toro: Riemann solvers,3rd, chap 17
   
    if(stp==1) {  
       if(sqrt(x0*x0+y0*y0)<= 0.4){
          r=1.0;
        } 
       else {
           r=0.125;
        }
   }
   else if(stp==2|| stp==3) {
       r=0.0;
   }
   else if(stp==4) {
      if(sqrt(x0*x0+y0*y0)<= 0.4){
         r=1.0/0.4;
       }
       else
       {r=0.1/0.4;} 
   }     
  */         
    return r;
} 

//void read_sol_data(string &infn);

void initialDataGrid()
{
  double xc, yc, dx, dy;
   
  Node *current;
  current = HeadListAllGrid;  
  OctCell *lsbl=NULL;    
  while(current != NULL)
  {		
    lsbl = current->cell;
    xc=lsbl->xc1;
    yc=lsbl->yc1;  
    dx=lsbl->dx;
    dy=lsbl->dy;
    lsbl->set_invM();     
    if(lsbl->flag % 2== 0 ) {   
      for(int imm=1; imm<=4; imm++)  {
        lsbl->dof0[imm-1][0]=lsbl->dof[imm-1][0]=u0xy(0.0, 0.0, xc,yc,dx,dy,imm);
        for(int i=1; i<=nDOF-1; i++) {
           lsbl->dof0[imm-1][i]=lsbl->dof[imm-1][i]=0.0;		    
         }	        
      }
    }   
    current = current->next;
 }
  //      string ifile("IN//in280.plt");
 //       read_sol_data(ifile);       

}




