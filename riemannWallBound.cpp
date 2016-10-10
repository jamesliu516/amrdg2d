
#include"non_uniform_grid.h"

const double gama=GAMMA;

const double  g1 = (gama - 1.0)/(2.0*gama);
const double  g2 = (gama + 1.0)/(2.0*gama);
const double  g3 = 2.0*gama/(gama - 1.0);
const double  g4 = 2.0/(gama - 1.0);
const double  g5 = 2.0/(gama + 1.0);
const double  g6 = (gama - 1.0)/(gama + 1.0);
const double  g7 = (gama - 1.0)/2.0;
const double  g8 = gama - 1.0;
//here n normal to wall(from fluid to wall)
void  solveERS(double dl,double ul, double pl, 
			  double uw, double *di,double *pi)
{
    double cM, A_M,B_M,aM;
    double du0;
    cM=ul-uw;
    A_M=g5/dl;
    B_M=pl*g6;

    du0=ul-uw;
    aM=sqrt(gama*pl/dl);

    if(du0<0.0){
        *pi=pl*pow(1.0+g7*cM/aM, g3);
        *di=dl*pow(*pi/pl, 1.0/gama);
    }
    else
    {
        *pi=pl+cM*0.5/A_M*(cM+sqrt(cM*cM+4.0*A_M*(B_M+pl)));
        *di=dl*(*pi/pl+g6)/(g6*(*pi)/pl+1.0);
    }
}

//hr 参考点到墙的距离，hg壁面到ghost点的距离
        
void reflectionRie(double hr,double hg, double vw,double U[4], PXYZ *nml, JBBL *jbugst)
{
    double tx,ty,nx,ny;
    nx=nml->x;
    ny=nml->y;
    tx=ny;
    ty=-nx;
    Vec2D ndir(nx,ny);
    Vec2D tdir(tx,ty);
    Vec2D ur(U[1],U[2]);
    double un1,ut1,un2,ut2;
    un1=ur.dot(ndir);
    ut1=ur.dot(tdir);
    un2=-(hg/hr)*un1;
    ut2=ut1;

    jbugst->u = (ny * ut2 - ty * un2)/(tx * ny - nx * ty);
    jbugst->v = (tx * un2 - nx * ut2)/(tx * ny - nx * ty);			
    double dw,pw;

    solveERS(U[0],-un1,U[3],vw,&dw,&pw);

 //   jbugst->p=U[3];
 //   jbugst->q=U[0];
    jbugst->p=pw+hg/hr*(pw-U[3]);
    jbugst->q=dw+hg/hr*(dw-U[0]);
  //  jbugst->p=pw;
   // jbugst->q=dw;
    if(jbugst->p<1.0e-13 || jbugst->q < 1.0e-13||
        pw<1.0e-13|| dw<1.0e-13){
        jbugst->p=U[3];
        jbugst->q=U[0];
    }
}


double getWallCurvatureRadius(const PXYZ *pt, PXYZ &ndir);
double getWallCurvatureRadius(const PXYZ *pt);

void reflectionRieQVLV(double hr,double hg, double vw,double U[4], PXYZ *nml, const PXYZ *wpt, JBBL *jbugst)
{
    double rd,tx,ty,nx,ny;
    nx=nml->x;
    ny=nml->y;
    tx=ny;
    ty=-nx;
    rd=getWallCurvatureRadius(wpt);
    Vec2D ndir(nx,ny);
    Vec2D tdir(tx,ty);
    Vec2D ur(U[1],U[2]);
    double un1,ut1,un2,ut2;
    un1=ur.dot(ndir);
    ut1=ur.dot(tdir);
    un2=-(hg/hr)*un1;
    ut2=ut1;

    jbugst->u = (ny * ut2 - ty * un2)/(tx * ny - nx * ty);
    jbugst->v = (tx * un2 - nx * ut2)/(tx * ny - nx * ty);			
    double dw,pw;

    solveERS(U[0],-un1,U[3],vw,&dw,&pw);

    jbugst->p=pw+hg*((pw-U[3])/hr-dw*ut2*ut2/rd);
    jbugst->q=dw+hg/hr*(dw-U[0]);
    if(jbugst->p<1.0e-13|| jbugst->q<1.0e-13
        || pw<1.0e-13||dw<1.0e-13){
        jbugst->p=U[3];
        jbugst->q=U[0];
    }
}



