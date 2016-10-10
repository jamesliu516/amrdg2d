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
	
void jbbl2shbl(double ujb[],double ush[]);
//(x,y) belong to [-1,1]x[-1,1] 在正规标准正方形元上定义
double u0xy(double xt, double yt, double xc, double yc, double dx, double dy, int stp)
{
    double r;

    double x0, y0, hxt0;
    x0=0.5*xt*dx+xc;
    y0=0.5*yt*dy+yc; 

    JBBL ujb[4];
    SHBL ush[4];
   
    JBBL jbu1;
    SHBL shu1;
    /*
//    unsteady shock vortex interaction 
	double ur,pr,qr,ul,pl,ql;
    State_l.q=1.0;
    State_l.p=1.0;
    State_l.u=1.1*sqrt(GAMMA);
    State_l.v=0.0;

	ql=jbu1.q=State_l.q;
	pl=jbu1.p=State_l.p;
	ul=jbu1.u=State_l.u;
    jbu1.v=0.0;

    jbblToShbl(&jbu1, &shu1);
    JBBL jbu0;
    SHBL shu0;
    right_shock_sol(ql,ul,pl,0.0,qr,ur,pr);
    State_r.q=qr;
    State_r.p=pr;
    jbu0.q=qr;
    State_r.u=jbu0.u=ur;
    State_r.v=jbu0.v=0.0;
    jbu0.p=pr;
    jbblToShbl(&jbu0, &shu0);
    double dtT,dU,tau,rxy,epsi=0.3,Radc=0.05,alph=0.204;
    double xo,yo;
    double gama=GAMMA;
    double U[4]={0.0}, U1[4]={0.0};
    xo=0.25,yo=0.5;
    rxy=(x0-xo)*(x0-xo)+(y0-yo)*(y0-yo);
    if(rxy>1e-6)rxy=sqrt(rxy);

    tau=rxy/Radc;
    dtT=-(gama-1.0)*epsi*epsi*exp(2.0*alph*(1.0-tau*tau))/(4.0*alph*gama);
    dU=epsi*tau*exp(alph*(1.0-tau*tau));
    //U[0]=pow(1.0+dtT,2.5);
    double Tl,Tr,Sl,Sr;
    Tl=pl/ql;
    Tr=pr/qr;
    Sl=pl/pow(ql,GAMMA);
    Sr=pr/pow(qr,GAMMA);
    U[0]=pow((Tl+dtT)/Sl,1.0/(GAMMA-1));
    U1[0]=pow((Tr+dtT)/Sr,1.0/(GAMMA-1));
    if(rxy>1e-6){
        U[1]=State_l.u+dU*(y0-yo)/rxy;
        U[2]=State_l.v-dU*(x0-xo)/rxy;
        U1[1]=State_r.u+dU*(y0-yo)/rxy;
        U1[2]=State_r.v-dU*(x0-xo)/rxy;
    }else{		
        U[1]=State_l.u;	
        U[2]=State_l.v;
        U1[1]=State_r.u;	
        U1[2]=State_r.v;
    }
  //  U[3]=pow(U[0],1.4);			
    U[3]=U[0]*(Tl+dtT);			
    U1[3]=U1[0]*(Tr+dtT);			
    double W[4]={0.0};
    double W1[4]={0.0};
    jbbl2shbl(U,W);
    jbbl2shbl(U1,W1);

    switch(stp){
    case 1:
        if(x0>0.5) r=W1[0];
        else r=W[0];
        break;
    case 2:
        if(x0>0.5)  r=W1[1];
        else  r=W[1];
        break;
    case 3:
        if(x0>0.5)  r=W1[2];
        else  r=W[2];
        break;
    case 4:
        if(x0>0.5)  r=W1[3];
        else  r=W[3];
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
//    unsteady shock triangle 
	double ur,pr,qr,ul,pl,ql;
    State_r.q=1.4;
    State_r.p=1.0;
    State_r.u=0.0;
    State_r.v=0.0;

	qr=jbu1.q=State_r.q;
	pr=jbu1.p=State_r.p;
	ur=jbu1.u=jbu1.v=0.0;

	jbblToShbl(&jbu1, &shu1);
    JBBL jbu0;
    SHBL shu0;
	left_shock_sol(qr,ur,pr,10.0,ql,ul,pl);
    State_l.q=ql;
    State_l.p=pl;
	jbu0.q=ql;
	State_l.u=jbu0.u=ul;
	State_l.v=jbu0.v=0.0;
	jbu0.p=pl;
	jbblToShbl(&jbu0, &shu0);
	switch(stp){
	case 1:
		if(x0<0.0) r=shu0.q;
		else r=shu1.q;
		break;
	case 2:
		if(x0<0.0)  r=shu0.qu;
		else  r=shu1.qu;
		break;
	case 3:
		if(x0<0.0)  r=shu0.qv;
		else  r=shu1.qv;
		break;
	case 4:
		if(x0<0.0)  r=shu0.te;
		else  r=shu1.te;
		break;
	}


   /*
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

double funcMult(double (* fun0)(double,double),
    double (*fun1)(double,double,double,double,double,double), double x,double y,double xc,double yc, double dx, double dy   )
{
    return (fun0(x,y)*fun1(x,y,xc,yc,dx,dy));
}

double ubase0(double x,double y,double xc,double yc, double dx, double dy)
{
    return 1.0;
}

double ubase1(double x,double y,double xc,double yc, double dx, double dy)
{
    return phix(x, xc, dx);
}

double ubase2(double x,double y,double xc,double yc, double dx, double dy)
{
    return psiy(y, yc, dy);
}

double ubase3(double x,double y,double xc,double yc, double dx, double dy)
{
    return phix_psiy(x,y, xc,yc,dx,dy);
}

double ubase4(double x,double y,double xc,double yc, double dx, double dy)
{
    return phix2m(x, xc, dx);
}

double ubase5(double x,double y,double xc,double yc, double dx, double dy)
{
    return psiy2m(y, yc, dy);
}

double u0xyGen(double x0, double y0, int stp);

double u0xyGen_q(double x,double y)
{
    return u0xyGen(x,y,1);
}
    
double u0xyGen_qu(double x,double y)
{
    return u0xyGen(x,y,2);
}

double u0xyGen_qv(double x,double y)
{
    return u0xyGen(x,y,3);
}

double u0xyGen_te(double x,double y)
{
    return u0xyGen(x,y,4);
}

double (*fun_u0Gen_ar[4])(double,double)={u0xyGen_q,u0xyGen_qu,
    u0xyGen_qv,u0xyGen_te};
//一般坐标
double u0xyGen(double x0, double y0, int stp)
{
    double r;

    double hxt0;

    JBBL ujb[4];
    SHBL ush[4];
   
    JBBL jbu1;
    SHBL shu1;
    /*
//    unsteady shock vortex interaction 
	double ur,pr,qr,ul,pl,ql;
    State_l.q=1.0;
    State_l.p=1.0;
    State_l.u=1.1*sqrt(GAMMA);
    State_l.v=0.0;

	ql=jbu1.q=State_l.q;
	pl=jbu1.p=State_l.p;
	ul=jbu1.u=State_l.u;
    jbu1.v=0.0;

    jbblToShbl(&jbu1, &shu1);
    JBBL jbu0;
    SHBL shu0;
    right_shock_sol(ql,ul,pl,0.0,qr,ur,pr);
    State_r.q=qr;
    State_r.p=pr;
    jbu0.q=qr;
    State_r.u=jbu0.u=ur;
    State_r.v=jbu0.v=0.0;
    jbu0.p=pr;
    jbblToShbl(&jbu0, &shu0);
    double dtT,dU,tau,rxy,epsi=0.3,Radc=0.05,alph=0.204;
    double xo,yo;
    double gama=GAMMA;
    double U[4]={0.0}, U1[4]={0.0};
    xo=0.25,yo=0.5;
    rxy=(x0-xo)*(x0-xo)+(y0-yo)*(y0-yo);
    if(rxy>1e-6)rxy=sqrt(rxy);

    tau=rxy/Radc;
    dtT=-(gama-1.0)*epsi*epsi*exp(2.0*alph*(1.0-tau*tau))/(4.0*alph*gama);
    dU=epsi*tau*exp(alph*(1.0-tau*tau));
    //U[0]=pow(1.0+dtT,2.5);
    double Tl,Tr,Sl,Sr;
    Tl=pl/ql;
    Tr=pr/qr;
    Sl=pl/pow(ql,GAMMA);
    Sr=pr/pow(qr,GAMMA);
    U[0]=pow((Tl+dtT)/Sl,1.0/(GAMMA-1));
    U1[0]=pow((Tr+dtT)/Sr,1.0/(GAMMA-1));
    if(rxy>1e-6){
        U[1]=State_l.u+dU*(y0-yo)/rxy;
        U[2]=State_l.v-dU*(x0-xo)/rxy;
        U1[1]=State_r.u+dU*(y0-yo)/rxy;
        U1[2]=State_r.v-dU*(x0-xo)/rxy;
    }else{		
        U[1]=State_l.u;	
        U[2]=State_l.v;
        U1[1]=State_r.u;	
        U1[2]=State_r.v;
    }
  //  U[3]=pow(U[0],1.4);			
    U[3]=U[0]*(Tl+dtT);			
    U1[3]=U1[0]*(Tr+dtT);			
    double W[4]={0.0};
    double W1[4]={0.0};
    jbbl2shbl(U,W);
    jbbl2shbl(U1,W1);

    switch(stp){
    case 1:
        if(x0>0.5) r=W1[0];
        else r=W[0];
        break;
    case 2:
        if(x0>0.5)  r=W1[1];
        else  r=W[1];
        break;
    case 3:
        if(x0>0.5)  r=W1[2];
        else  r=W[2];
        break;
    case 4:
        if(x0>0.5)  r=W1[3];
        else  r=W[3];
        break;
    }
*/
    //1233487494w4852545234545446jje5h45thgehrrhhy
//    unsteady shock triangle
	double ur,pr,qr,ul,pl,ql;
    double ma_free=10.0;
    State_r.q=1.4;
    State_r.p=1.0;
    State_r.u=0.0;
    State_r.v=0.0;

	qr=jbu1.q=State_r.q;
	pr=jbu1.p=State_r.p;
	ur=jbu1.u=jbu1.v=0.0;

	jbblToShbl(&jbu1, &shu1);
    JBBL jbu0;
    SHBL shu0;
	left_shock_sol(qr,ur,pr,ma_free,ql,ul,pl);
    State_l.q=ql;
    State_l.p=pl;
	jbu0.q=ql;
	State_l.u=jbu0.u=ul;
	State_l.v=jbu0.v=0.0;
	jbu0.p=pl;
	jbblToShbl(&jbu0, &shu0);
	switch(stp){
	case 1:
		if(x0<-0.2) r=shu0.q;
		else r=shu1.q;
		break;
	case 2:
		if(x0<-0.2)  r=shu0.qu;
		else  r=shu1.qu;
		break;
	case 3:
		if(x0<-0.2)  r=shu0.qv;
		else  r=shu1.qv;
		break;
	case 4:
		if(x0<-0.2)  r=shu0.te;
		else  r=shu1.te;
		break;
	}
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
   /*
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

double (*fun_ar[6])(double ,double ,double ,double , double , double )={ubase0,ubase1,ubase2,ubase3,ubase4,ubase5};
double funcMult(double (* fun0)(double,double),
    double (*fun1)(double,double,double,double,double,double), double x,double y,double xc,double yc, double dx, double dy);
void initialDataGrid()
{
    double xc, yc, dx, dy;

    double gs4pts[4]={-0.861136311594053,-0.339981043584856,0.339981043584856,0.861136311594053};
    double wei4pts[4]={0.347854845137454, 0.652145154862546, 0.652145154862546,0.347854845137454};
    Node *current;
    current = HeadListAllGrid;  
    OctCell *lsbl=NULL;    
    double jfSum=0.0;
    while(current != NULL)
    {		
        lsbl = current->cell;
        xc=lsbl->xc1;
        yc=lsbl->yc1;  
        dx=lsbl->dx;
        dy=lsbl->dy;
        lsbl->set_invM();     
 //       if(lsbl->flag % 2== 0 ) {   
        if(1 ) {   
            for(int imm=1; imm<=4; imm++)  {
                for(int i=0; i<=nDOF-1; i++) {
                    jfSum=0.0;
                    for(int m0=0;m0<4;++m0){
                        for(int m1=0;m1<4;++m1){
                            jfSum+=wei4pts[m0]*wei4pts[m1]*
                                funcMult(fun_u0Gen_ar[imm-1],fun_ar[i],
                                    xc+0.5*dx*gs4pts[m1], yc+0.5*dy*gs4pts[m0],
                                    xc,yc,dx,dy);
                        }
                    }
                    jfSum*=0.5*dx*dy;

                    lsbl->dof0[imm-1][i]=lsbl->dof[imm-1][i]=
                        lsbl->invM[i]*jfSum;		    
                }	        
            }
            //----
            for(int imm=1; imm<=4; imm++)  {
                lsbl->dof0[imm-1][0]=lsbl->dof[imm-1][0]=u0xy(0.0, 0.0, xc,yc,dx,dy,imm);
                for(int i=1; i<=nDOF-1; i++) {
                    lsbl->dof0[imm-1][i]=lsbl->dof[imm-1][i]=0.0;		    
                }	        
            }
            //=====
        }   
        current = current->next;
    }
    //      string ifile("IN//in280.plt");
    //       read_sol_data(ifile);       

    }
