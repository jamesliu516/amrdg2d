/*********************************************************888********
//如果尖后缘在一个网格中,退化的网格那这断尖后缘被截断,网格当着计算单元
//多值点注意这个尖后缘点还有前缘点
//需要找到那个前缘点,尖后援所在的网格,那个网格特殊处理.有一个face的ghost
//cell需要平均上下面ghost cell的结果
//这个网格可能会标记为界面网格,或者就有可能就标记为流体网格
//现在尖后缘被截断,那这个网格就是计算单元
//如果面两端有一个是需要处理的固体内部点,那所对应的ghost点与对称点的对 
//称线以离按面扫描的面中点最近的那个翼形表面点的两端的两个折线断的圆
//作为径向对称线
//边界法线的处理这里借鉴dg的一篇文章关于高阶格式低阶边界近似的处理方法
//关于多值点最后处理通过适当的标记指比如标志为4来特殊处理
 *********************************************************************/
#include"non_uniform_grid.h"
#include"findneighbor.h"
#include<set>
#include<string>
#include<cstddef>
#include<utility>
#include<algorithm>

extern OctCell *bodygrid;
extern BndryNode *HeadFaceCell;
extern vector<OctCell*> ExtWallBndry;

extern Node *HeadListAllGrid; //the head of all cell 
void shblToJbbl(SHBL *ush, JBBL *jbu);
void jbblToShbl(JBBL *jbu, SHBL *ush);

//double geometry_fun(double,double);
//double distance2p(const PXYZ &p1, const PXYZ &p2);

int getIndexAtWallPoint(const PXYZ *gst, int &nb);
//extern PXYZ Circle[NPOINTCIRCLE+1];

//WALLP Company[NPOINTCIRCLE+1]; //壁面输出量
WALLP MultiCompany[N_BODY+1][500]; //壁面输出量
extern PXYZ MultiCircle[N_BODY+1][500]; //points at solid wall
extern int NWallPts[N_BODY+3];

extern PXYZ ext_wall[NPT_EXT_WALL+1]; //points at external computational domain
//const PXYZ *const wallpoint = Circle;

const int NumberAirfoilPoint=NPOINTCIRCLE;

//const int HalfNumberAirfoilPoint=NumberAirfoilPoint / 2;

//得到AB连线为对称线的对称点
void getsmtrctAB(const PXYZ *gst, const PXYZ *pA, const PXYZ *pB,
    PXYZ *smtrctP, PXYZ *jiaodian, PXYZ *normal)
{
    double ls1;
    double ls2;
    double b1, b2, xjd, yjd, jl;

    if(fabs(pB->y-pA->y) < 1e-10)
    {
        smtrctP->x = gst->x;
        smtrctP->y = 2.0 * pA->y - gst->y;

        normal->x = 0.0;
        normal->y = (smtrctP->y - gst->y)/fabs(smtrctP->y - gst->y);

        jiaodian->x = gst->x;
        jiaodian->y = pA->y;
    }
    else if(fabs(pA->x - pB->x) < 1e-10)
    {
        smtrctP->x = 2 * pA->x - gst->x;
        smtrctP->y = gst->y;

        jiaodian->x = pA->x;
        xjd=jiaodian->x;

        jiaodian->y = gst->y;
        yjd=jiaodian->y;

        jl = sqrt((gst->x - xjd) * (gst->x - xjd) + (gst->y - yjd) * (gst->y - yjd));

        normal->x = (xjd - gst->x)/jl;
        normal->y = (yjd - gst->y)/jl;
    }
    else
    {
        ls1 = (pA->y-pB->y)/(pA->x-pB->x);
        ls2 = - 1.0/ls1;
        b1 = pA->y - ls1 * pA->x;
        b2 = gst->y - ls2 * gst->x;

        xjd = (b2 - b1)/(ls1-ls2);
        yjd = ls1 * xjd + b1;

        jiaodian->x = xjd;
        jiaodian->y = yjd;

        smtrctP->x = 2 * xjd - gst->x;
        smtrctP->y = 2 * yjd - gst->y;	

        jl = sqrt((gst->x-xjd) * (gst->x-xjd) + (gst->y-yjd) * (gst->y-yjd));

        normal->x = (xjd - gst->x)/jl;
        normal->y = (yjd - gst->y)/jl;
    }
}

void getsmtrctAB(double hd, const PXYZ *gst, const PXYZ *pA, const PXYZ *pB,
    PXYZ *smtrctP, PXYZ *jiaodian, PXYZ *normal, double *dis)
{
    double ls1;
    double ls2;
    double b1, b2, xjd, yjd, jl;

    if(fabs(pB->y-pA->y) < 1e-10)
    {
//        smtrctP->x = gst->x;
  //      smtrctP->y = 2.0 * pA->y - gst->y;

        normal->x = 0.0;
        normal->y = (smtrctP->y - gst->y)/fabs(smtrctP->y - gst->y);
        xjd=jiaodian->x = gst->x;
        yjd=jiaodian->y = pA->y;
        *dis = sqrt((gst->x-xjd) * (gst->x-xjd) + (gst->y-yjd) * (gst->y-yjd));
        smtrctP->x=xjd+hd*normal->x;
        smtrctP->y=yjd+hd*normal->y;
    }
    else if(fabs(pA->x - pB->x) < 1e-10)
    {
       // smtrctP->x = 2 * pA->x - gst->x;
        //smtrctP->y = gst->y;

        jiaodian->x = pA->x;
        xjd=jiaodian->x;

        jiaodian->y = gst->y;
        yjd=jiaodian->y;

        *dis=jl = sqrt((gst->x - xjd) * (gst->x - xjd) + (gst->y - yjd) * (gst->y - yjd));

        normal->x = (xjd - gst->x)/jl;
        normal->y = (yjd - gst->y)/jl;
        smtrctP->x=xjd+hd*normal->x;
        smtrctP->y=yjd+hd*normal->y;
    }
    else
    {
        ls1 = (pA->y-pB->y)/(pA->x-pB->x);
        ls2 = - 1.0/ls1;
        b1 = pA->y - ls1 * pA->x;
        b2 = gst->y - ls2 * gst->x;

        xjd = (b2 - b1)/(ls1-ls2);
        yjd = ls1 * xjd + b1;

        jiaodian->x = xjd;
        jiaodian->y = yjd;

//        smtrctP->x = 2 * xjd - gst->x;
  //      smtrctP->y = 2 * yjd - gst->y;	

        *dis=jl = sqrt((gst->x-xjd) * (gst->x-xjd) + (gst->y-yjd) * (gst->y-yjd));

        normal->x = (xjd - gst->x)/jl;
        normal->y = (yjd - gst->y)/jl;
        smtrctP->x=xjd+hd*normal->x;
        smtrctP->y=yjd+hd*normal->y;
    }
}

/*先要找到对称点处的jbu然后根据此函数就可以得到ghost点值, nxyz为法向向量*/
void reflection(const JBBL *jbu, const PXYZ *nxyz, JBBL *jbugst)
{
    jbugst->q = jbu->q;
    jbugst->p = jbu->p;

    jbugst->u = (- nxyz->x * nxyz->x + nxyz->y * nxyz->y) * jbu->u - 2 * nxyz->x * nxyz->y * jbu->v; 
    jbugst->v = - 2 * nxyz->x * nxyz->y * jbu->u + (nxyz->x * nxyz->x - nxyz->y * nxyz->y)* jbu->v; 
}
//hr 参考点到墙的距离，hg壁面到ghost点的距离
void reflection(double hr, double hg,const JBBL *jbu, const PXYZ *nxyz, JBBL *jbugst)
{
    double tx,ty,nx,ny;
    nx=nxyz->x;
    ny=nxyz->y;
    tx=ny;
    ty=-nx;
    Vec2D ndir(nx,ny);
    Vec2D tdir(tx,ty);
    Vec2D ur(jbu->u,jbu->v);
    double un1,ut1,un2,ut2;
    un1=ur.dot(ndir);
    ut1=ur.dot(tdir);
    un2=-(hg/hr)*un1;
    ut2=ut1;

    jbugst->q = jbu->q;
    jbugst->p = jbu->p;
    jbugst->u = (ny * ut2 - ty * un2)/(tx * ny - nx * ty);
    jbugst->v = (tx * un2 - nx * ut2)/(tx * ny - nx * ty);			
}

double getWallCurvatureRadius(const PXYZ *pt);
double getWallCurvatureRadius(const PXYZ *pt, PXYZ &ndir);
//使用非镜像点使用参考点的处理
/*
void reflectionForMoving6(double hr,double hg, const JBBL *jbu,				
    const PXYZ *wall, const PXYZ *gstc, const PXYZ *nxyz, 
    double vnsolid, JBBL *jbugst)
{
    double tx,ty;
    double jl, un1,ut1, un2,ut2, nx, ny;
    double sign, rd;
#ifdef NORMAL_DIR_MOD
    PXYZ ndir;
   // jl = sqrt((wall->x - gstc->x) * (wall->x - gstc->x)
     //   +(wall->y - gstc->y) * (wall->y - gstc->y));
    rd=getWallCurvatureRadius(wall, ndir);
    tx=ndir.y;
    ty=-ndir.x;
    un1=jbu->u * ndir.x + jbu->v * ndir.y;
    ut1=jbu->u * tx + jbu->v * ty;			
    if(ut1>=0.0) sign = 1.0;
    else sign = -1.0;
    jbugst->p = jbu->p - jbu->q * ut1 * ut1 * (hr+hg)/rd; //曲率半径符号问题,这里只考虑了曲率中心在固体中情形

    jbugst->q = jbu->q * pow(jbugst->p/jbu->p, 1.0/GAMMA);	
    ut2=sign * sqrt(ut1*ut1 + (2*GAMMA/(GAMMA-1)) * (jbu->p/jbu->q - jbugst->p/jbugst->q));

    un2= - hg/hr*un1;

    nx=ndir.x;
    ny=ndir.y;
    jbugst->u = (ny * ut2 - ty * un2)/(tx * ny - nx * ty);
    jbugst->v = (tx * un2 - nx * ut2)/(tx * ny - nx * ty);			
#else	

    //	jbugst->q = jbu->q;
    //	jbugst->p = jbu->p;
    //	hc=h/pow(2.0, Nr);
   // jl = sqrt((wall->x - gstc->x) * (wall->x - gstc->x)
     //   +(wall->y - gstc->y) * (wall->y - gstc->y));
    nx=nxyz->x;
    ny=nxyz->y;

    rd=getWallCurvatureRadius(wall);

    //   cout<<"cur R=" << rd << endl;

    tx=ny;
    ty=-nx;
    un1=jbu->u * nx + jbu->v * ny;
    ut1=jbu->u * tx + jbu->v * ty;
    if(ut1>=0.0) sign = 1.0;
    else sign = -1.0;

    //	jbugst->p = 0.5*(jbu->p + pw + jl * (pw - pwh)/hc);
    jbugst->p = jbu->p - jbu->q * ut1 * ut1 * (hr+hg)/rd; //曲率半径符合问题,这里只考虑了曲率中心在固体中情形

    jbugst->q = jbu->q * pow(jbugst->p/jbu->p, 1.0/GAMMA);
    //	jbugst->q = qw1 + jl *(qw1-qwh1)/hc;
    //-------------------------------------------------------------
    //	jbugst->u = (- nxyz->x * nxyz->x + nxyz->y * nxyz->y) * jbu->u
     //  - 2 * nxyz->x * nxyz->y * jbu->v + 2 * nxyz->x * vnsolid;
    //    jbugst->v = - 2 * nxyz->x * nxyz->y * jbu->u + (nxyz->x * nxyz->x - nxyz->y * nxyz->y)* jbu->v 
//        + 2 * nxyz->y * vnsolid; 

    ut2=sign * sqrt(ut1*ut1 + (2*GAMMA/(GAMMA-1)) * (jbu->p/jbu->q - jbugst->p/jbugst->q));

    un2= - hg/hr*un1;

    jbugst->u = (ny * ut2 - ty * un2)/(tx * ny - nx * ty);
    jbugst->v = (tx * un2 - nx * ut2)/(tx * ny - nx * ty);
#endif
}
*/


void reflectionForMoving6(double hc, const JBBL *jbu,				
    const PXYZ *wall, const PXYZ *gstc, const PXYZ *nxyz, 
    double vnsolid, JBBL *jbugst)
{
    double tx,ty;
    double jl, un1,ut1, un2,ut2, nx, ny;
    double sign, rd;
    double dil,dir,ui,prei;
#ifdef NORMAL_DIR_MOD
    PXYZ ndir;
    jl = sqrt((wall->x - gstc->x) * (wall->x - gstc->x)
        +(wall->y - gstc->y) * (wall->y - gstc->y));
    rd=getWallCurvatureRadius(wall, ndir);
    if(rd>1.0e10){
        ndir.x=nxyz->x;
        ndir.y=nxyz->y;
    }
    tx=ndir.y;
    ty=-ndir.x;
    un1=jbu->u * ndir.x + jbu->v * ndir.y;
    ut1=jbu->u * tx + jbu->v * ty;			
    if(ut1>=0.0) sign = 1.0;
    else sign = -1.0;
    jbugst->p = jbu->p - jbu->q * ut1 * ut1 * 2 *jl/rd; //曲率半径符号问题,这里只考虑了曲率中心在固体中情形

    jbugst->q = jbu->q * pow(jbugst->p/jbu->p, 1.0/GAMMA);	
   // ut2=sign * sqrt(ut1*ut1 + (2*GAMMA/(GAMMA-1)) * (jbu->p/jbu->q - jbugst->p/jbugst->q));
    ut2=ut1;

    un2= - un1;
/*    if(b_ARS){
        solveARS(jbugst->q,jbu->q,1.4,1.4,0.0,0.0,un2,un1,jbugst->p,jbu->p,
            &dil,&dir,&ui,&prei);
        //cout<<fabs(dil-dir)<<" "<<fabs(prei-jbugst->p)<<endl;
        jbugst->p=prei;
     //   jbugst->q=0.5*(dil+dir);
        jbugst->q = jbu->q * pow(jbugst->p/jbu->p, 1.0/GAMMA);
    }*/
    nx=ndir.x;
    ny=ndir.y;
    jbugst->u = (ny * ut2 - ty * un2)/(tx * ny - nx * ty);
    jbugst->v = (tx * un2 - nx * ut2)/(tx * ny - nx * ty);			
#else	

    //	jbugst->q = jbu->q;
    //	jbugst->p = jbu->p;
    //	hc=h/pow(2.0, Nr);
    jl = sqrt((wall->x - gstc->x) * (wall->x - gstc->x)
        +(wall->y - gstc->y) * (wall->y - gstc->y));
    nx=nxyz->x;
    ny=nxyz->y;

    rd=getWallCurvatureRadius(wall);

    //   cout<<"cur R=" << rd << endl;

    tx=ny;
    ty=-nx;
    un1=jbu->u * nx + jbu->v * ny;
    ut1=jbu->u * tx + jbu->v * ty;
    if(ut1>=0.0) sign = 1.0;
    else sign = -1.0;

    //	jbugst->p = 0.5*(jbu->p + pw + jl * (pw - pwh)/hc);
    jbugst->p = jbu->p - jbu->q * ut1 * ut1 * 2 *jl/rd; //曲率半径符合问题,这里只考虑了曲率中心在固体中情形
 if(jbu->p<1e-10)jbugst->p = jbu->p ;
    jbugst->q = jbu->q * pow(jbugst->p/jbu->p, 1.0/GAMMA);
    //	jbugst->q = qw1 + jl *(qw1-qwh1)/hc;
    //-------------------------------------------------------------
    /*	jbugst->u = (- nxyz->x * nxyz->x + nxyz->y * nxyz->y) * jbu->u
        - 2 * nxyz->x * nxyz->y * jbu->v + 2 * nxyz->x * vnsolid; 

        jbugst->v = - 2 * nxyz->x * nxyz->y * jbu->u + (nxyz->x * nxyz->x - nxyz->y * nxyz->y)* jbu->v 
        + 2 * nxyz->y * vnsolid; */

//    ut2=sign * sqrt(ut1*ut1 + (2*GAMMA/(GAMMA-1)) * (jbu->p/jbu->q - jbugst->p/jbugst->q));
    ut2=ut1;
    un2= - un1;
    /*if(b_ARS){
        solveARS(jbugst->q,jbu->q,1.4,1.4,0.0,0.0,un2,un1,jbugst->p,jbu->p,
            &dil,&dir,&ui,&prei);
        //cout<<fabs(dil-dir)<<" "<<fabs(prei-jbugst->p)<<endl;
        jbugst->p=prei;
    //    jbugst->q=0.5*(dil+dir);
    jbugst->q = jbu->q * pow(jbugst->p/jbu->p, 1.0/GAMMA);
    }*/
    //void  solveARS(double dl,double dr,double gl,double gr,double bl,
    //			   double br,double ul,double ur,double pl,double pr,
    //			   double *dil,double *dir,double *ui,double *pi);

    jbugst->u = (ny * ut2 - ty * un2)/(tx * ny - nx * ty);
    jbugst->v = (tx * un2 - nx * ut2)/(tx * ny - nx * ty);

    //	jbugst->u = (- nxyz->x * nxyz->x + nxyz->y * nxyz->y) * jbu->u - 2 * nxyz->x * nxyz->y * jbu->v; 
    //	jbugst->v = - 2 * nxyz->x * nxyz->y * jbu->u + (nxyz->x * nxyz->x - nxyz->y * nxyz->y)* jbu->v; 

#endif
}

/*
   void getsmtrctAB(const PXYZ *gst, const PXYZ *pA, const PXYZ *pB,
   PXYZ *smtrctP, PXYZ *jiaodian, PXYZ *normal)*/
//得到AB连线为对称线的对称点
void getInterAB(const PXYZ *gst, const PXYZ *pA, const PXYZ *pB, PXYZ *jiaodian)
{
    double ls1;
    double ls2;
    double b1, b2, xjd, yjd, jl;

    if(fabs(pB->y-pA->y) < 1e-10)
    {
        jiaodian->x = gst->x;
        jiaodian->y = pA->y;
    }
    else if(fabs(pA->x - pB->x) < 1e-10)
    {
        jiaodian->x = pA->x;
        jiaodian->y = gst->y;
    }
    else
    {
        ls1 = (pA->y-pB->y)/(pA->x-pB->x);
        ls2 = - 1.0/ls1;
        b1 = pA->y - ls1 * pA->x;
        b2 = gst->y - ls2 * gst->x;

        xjd = (b2 - b1)/(ls1-ls2);
        yjd = ls1 * xjd + b1;

        jiaodian->x = xjd;
        jiaodian->y = yjd;
    }
}


//comp_c:表示网格中心在固体内部，且至少一个邻居网格的中心落在流场里面的
//引入用来解决TODO的问题
//int getIndexAtWallPoint(const PXYZ *gst,OctCell *comp_c=NULL)
/*
int getIndexAtWallPoint(const PXYZ *gst)
{
    int ils, ia,ib,i;
    double dis,dis2;
    PXYZ pa,pi1,pi2,pb1,pb2;

    ia=0;
    //	ib=NumberAirfoilPoint;

    dis=(gst->x-wallpoint[0].x)*(gst->x-wallpoint[0].x)
        +(gst->y-wallpoint[0].y)*(gst->y-wallpoint[0].y);

    for(i=1;i<=NumberAirfoilPoint-1;++i)
    {
        dis2 = (gst->x-wallpoint[i].x)*(gst->x-wallpoint[i].x)
            +(gst->y-wallpoint[i].y)*(gst->y-wallpoint[i].y);
        if(dis2<dis){
            dis=dis2;
            ia=i;
        }
    }
    if(ia==0)
    {
        pa.x = wallpoint[ia].x;
        pa.y = wallpoint[ia].y;
        pb1.x = wallpoint[ia+1].x; 
        pb1.y = wallpoint[ia+1].y;		

        pb2.x = wallpoint[NumberAirfoilPoint-1].x;
        pb2.y = wallpoint[NumberAirfoilPoint-1].y;
        getInterAB(gst, &pa, &pb1, &pi1);
        getInterAB(gst, &pa, &pb2, &pi2);

        //	if((pi1.x-pa.x)*(pb1.x-pa.x) + (pi1.y-pa.y)*(pb1.y-pa.y) >= 0.0) ia=ia; //2011-12-23修改,这里错了
        //	else if((pi2.x-pa.x)*(pb2.x-pa.x) + (pi2.y-pa.y)*(pb2.y-pa.y) >= 0.0) ia=NumberAirfoilPoint-1;
        if((pa.x-pi1.x)*(pb1.x-pi1.x) + (pa.y-pi1.y)*(pb1.y-pi1.y) <= 0.0) ia=ia; //2011-12-23修改,这里错了
        else if((pa.x-pi2.x)*(pb2.x-pi2.x) + (pa.y-pi2.y)*(pb2.y-pi2.y) <= 0.0) ia=NumberAirfoilPoint-1;
        else ia=ia;	
    }
    else
    {
        pa.x = wallpoint[ia].x;
        pa.y = wallpoint[ia].y;
        pb1.x = wallpoint[ia+1].x; 
        pb1.y = wallpoint[ia+1].y;		

        pb2.x = wallpoint[ia-1].x;
        pb2.y = wallpoint[ia-1].y;
        getInterAB(gst, &pa, &pb1, &pi1);
        getInterAB(gst, &pa, &pb2, &pi2);

        //	if((pi1.x-pa.x)*(pb1.x-pa.x) + (pi1.y-pa.y)*(pb1.y-pa.y) >= 0.0) ia=ia;
        //	else if((pi2.x-pa.x)*(pb2.x-pa.x) + (pi2.y-pa.y)*(pb2.y-pa.y) >= 0.0) ia=ia-1;

        if((pa.x-pi1.x)*(pb1.x-pi1.x) + (pa.y-pi1.y)*(pb1.y-pi1.y) <= 0.0) ia=ia; //2011-12-23修改,这里错了
        else if((pa.x-pi2.x)*(pb2.x-pi2.x) + (pa.y-pi2.y)*(pb2.y-pi2.y) <= 0.0) ia=ia-1;
        else ia=ia;		
    }
    return ia;
}
*/
//nb 指属于哪个物体
int getIndexAtWallPoint(const PXYZ *gst, int &nb)
{
    int ils, ia,ib,i;
    double dis,dis2;
    PXYZ pa,pi1,pi2,pb1,pb2;
    int wy;

    //	ib=NumberAirfoilPoint;
    dis=1.0e10;
    //  ia=0;
    for(wy=0; wy<N_BODY;++wy){
        for(i=0;i<=NWallPts[wy]-1;++i)
        {
            dis2 = (gst->x-MultiCircle[wy][i].x)*(gst->x-MultiCircle[wy][i].x)
                +(gst->y-MultiCircle[wy][i].y)*(gst->y-MultiCircle[wy][i].y);
            if(dis2<dis){
                dis=dis2;
                ia=i;
                nb=wy;
            }
        }
    }

    if(ia==0)
    {
        pa.x = MultiCircle[nb][ia].x;
        pa.y = MultiCircle[nb][ia].y;
        pb1.x = MultiCircle[nb][ia+1].x; 
        pb1.y = MultiCircle[nb][ia+1].y;		

        pb2.x = MultiCircle[nb][NWallPts[nb]-1].x;
        pb2.y = MultiCircle[nb][NWallPts[nb]-1].y;
        getInterAB(gst, &pa, &pb1, &pi1);
        getInterAB(gst, &pa, &pb2, &pi2);

        //	if((pi1.x-pa.x)*(pb1.x-pa.x) + (pi1.y-pa.y)*(pb1.y-pa.y) >= 0.0) ia=ia; //2011-12-23修改,这里错了
        //	else if((pi2.x-pa.x)*(pb2.x-pa.x) + (pi2.y-pa.y)*(pb2.y-pa.y) >= 0.0) ia=NumberAirfoilPoint-1;
        if((pa.x-pi1.x)*(pb1.x-pi1.x) + (pa.y-pi1.y)*(pb1.y-pi1.y) <= 0.0) ia=ia; //2011-12-23修改,这里错了
        else if((pa.x-pi2.x)*(pb2.x-pi2.x) + (pa.y-pi2.y)*(pb2.y-pi2.y) <= 0.0) ia=NWallPts[nb]-1;
        else ia=ia;	
    }
    else
    {
        pa.x = MultiCircle[nb][ia].x;
        pa.y = MultiCircle[nb][ia].y;
        pb1.x = MultiCircle[nb][ia+1].x; 
        pb1.y = MultiCircle[nb][ia+1].y;		

        pb2.x = MultiCircle[nb][ia-1].x;
        pb2.y = MultiCircle[nb][ia-1].y;
        getInterAB(gst, &pa, &pb1, &pi1);
        getInterAB(gst, &pa, &pb2, &pi2);

        //	if((pi1.x-pa.x)*(pb1.x-pa.x) + (pi1.y-pa.y)*(pb1.y-pa.y) >= 0.0) ia=ia;
        //	else if((pi2.x-pa.x)*(pb2.x-pa.x) + (pi2.y-pa.y)*(pb2.y-pa.y) >= 0.0) ia=ia-1;

        if((pa.x-pi1.x)*(pb1.x-pi1.x) + (pa.y-pi1.y)*(pb1.y-pi1.y) <= 0.0) ia=ia; //2011-12-23修改,这里错了
        else if((pa.x-pi2.x)*(pb2.x-pi2.x) + (pa.y-pi2.y)*(pb2.y-pi2.y) <= 0.0) ia=ia-1;
        else ia=ia;		
    }
    return ia;
}

int getIndexAtWallPoint(const PXYZ *gst, int NumberPoint,PXYZ wallpnt[]   )
{
    int ils, ia,ib,i;
    double dis,dis2;
    PXYZ pa,pi1,pi2,pb1,pb2;

    ia=0;
    //	ib=NumberAirfoilPoint;

    dis=(gst->x-wallpnt[0].x)*(gst->x-wallpnt[0].x)
        +(gst->y-wallpnt[0].y)*(gst->y-wallpnt[0].y);

    for(i=1;i<=NumberPoint-1;++i)
    {
        dis2 = (gst->x-wallpnt[i].x)*(gst->x-wallpnt[i].x)
            +(gst->y-wallpnt[i].y)*(gst->y-wallpnt[i].y);
        if(dis2<dis){
            dis=dis2;
            ia=i;
        }
    }
    if(ia==0)
    {
        pa.x = wallpnt[ia].x;
        pa.y = wallpnt[ia].y;
        pb1.x = wallpnt[ia+1].x; 
        pb1.y = wallpnt[ia+1].y;		

        pb2.x = wallpnt[NumberPoint-1].x;
        pb2.y = wallpnt[NumberPoint-1].y;
        getInterAB(gst, &pa, &pb1, &pi1);
        getInterAB(gst, &pa, &pb2, &pi2);

        //	if((pi1.x-pa.x)*(pb1.x-pa.x) + (pi1.y-pa.y)*(pb1.y-pa.y) >= 0.0) ia=ia; //2011-12-23修改,这里错了
        //	else if((pi2.x-pa.x)*(pb2.x-pa.x) + (pi2.y-pa.y)*(pb2.y-pa.y) >= 0.0) ia=NumberAirfoilPoint-1;
        if((pa.x-pi1.x)*(pb1.x-pi1.x) + (pa.y-pi1.y)*(pb1.y-pi1.y) <= 0.0) ia=ia; //2011-12-23修改,这里错了
        else if((pa.x-pi2.x)*(pb2.x-pi2.x) + (pa.y-pi2.y)*(pb2.y-pi2.y) <= 0.0) ia=NumberPoint-1;
        else ia=ia;	
    }
    else
    {
        pa.x = wallpnt[ia].x;
        pa.y = wallpnt[ia].y;
        pb1.x = wallpnt[ia+1].x; 
        pb1.y = wallpnt[ia+1].y;		

        pb2.x = wallpnt[ia-1].x;
        pb2.y = wallpnt[ia-1].y;
        getInterAB(gst, &pa, &pb1, &pi1);
        getInterAB(gst, &pa, &pb2, &pi2);

        //	if((pi1.x-pa.x)*(pb1.x-pa.x) + (pi1.y-pa.y)*(pb1.y-pa.y) >= 0.0) ia=ia;
        //	else if((pi2.x-pa.x)*(pb2.x-pa.x) + (pi2.y-pa.y)*(pb2.y-pa.y) >= 0.0) ia=ia-1;

        if((pa.x-pi1.x)*(pb1.x-pi1.x) + (pa.y-pi1.y)*(pb1.y-pi1.y) <= 0.0) ia=ia; //2011-12-23修改,这里错了
        else if((pa.x-pi2.x)*(pb2.x-pi2.x) + (pa.y-pi2.y)*(pb2.y-pi2.y) <= 0.0) ia=ia-1;
        else ia=ia;		
    }
    return ia;
}

double getRadius(const PXYZ *pt1, const PXYZ *pt2,const PXYZ *pt3)
{
    double rd,dls1,dls2,dls3,dls4;
    double a1,b1,a2,b2;
    dls1=(pt1->x-pt2->x)*(pt1->x-pt2->x)+(pt1->y-pt2->y)*(pt1->y-pt2->y);
    dls2=(pt1->x-pt3->x)*(pt1->x-pt3->x)+(pt1->y-pt3->y)*(pt1->y-pt3->y);
    dls3=(pt2->x-pt3->x)*(pt2->x-pt3->x)+(pt2->y-pt3->y)*(pt2->y-pt3->y);
    dls4=2.0*(pt3->x * (pt2->y - pt1->y) + pt2->x * (pt1->y - pt3->y) 
        + pt1->x * (pt3->y-pt2->y));
    a1=pt2->x-pt1->x;
    b1=pt2->y-pt1->y;

    a2=pt3->x-pt1->x;
    b2=pt3->y-pt1->y;
    rd = sqrt(dls1*dls2*dls3)/(fabs(dls4)+1.0e-18);
    if(fabs(a1*b2-b1*a2)<1.0e-12) rd=1.0e15;
    return rd;
}

double getWallCurvatureRadius(const PXYZ *pt)
{
    int ils1;
    double rd1=0.0,rd2=0.0;
    PXYZ pp1,pp2,pp3;
    int nb;
    ils1=getIndexAtWallPoint(pt,nb);
    if(ils1==NWallPts[nb]-1)
    {
        pp1.x = MultiCircle[nb][0].x;
        pp1.y = MultiCircle[nb][0].y;
        pp2.x = MultiCircle[nb][ils1].x;
        pp2.y = MultiCircle[nb][ils1].y;
        pp3.x = MultiCircle[nb][ils1-1].x;
        pp3.y = MultiCircle[nb][ils1-1].y;

        rd1=getRadius(&pp1,&pp2,&pp3);

        return rd1;
    }
    else if(ils1==NWallPts[nb]-2)
    {
        pp1.x = MultiCircle[nb][0].x;
        pp1.y = MultiCircle[nb][0].y;
        pp2.x = MultiCircle[nb][ils1].x;
        pp2.y = MultiCircle[nb][ils1].y;
        pp3.x = MultiCircle[nb][ils1+1].x;
        pp3.y = MultiCircle[nb][ils1+1].y;

        rd1=getRadius(&pp1,&pp2,&pp3);

        pp1.x = MultiCircle[nb][ils1+1].x;
        pp1.y = MultiCircle[nb][ils1+1].y;
        pp2.x = MultiCircle[nb][ils1].x;
        pp2.y = MultiCircle[nb][ils1].y;
        pp3.x = MultiCircle[nb][ils1-1].x;
        pp3.y = MultiCircle[nb][ils1-1].y;

        rd2=getRadius(&pp1,&pp2,&pp3);

        return (0.5*(rd1+rd2));
    }
    else if(ils1==0)
    {
        pp1.x = MultiCircle[nb][0].x;
        pp1.y = MultiCircle[nb][0].y;
        pp2.x = MultiCircle[nb][1].x;
        pp2.y = MultiCircle[nb][1].y;
        pp3.x = MultiCircle[nb][2].x;
        pp3.y = MultiCircle[nb][2].y;
        rd1=getRadius(&pp1,&pp2,&pp3);

        return rd1;
    }
    else
    {
        pp1.x = MultiCircle[nb][ils1+2].x;
        pp1.y = MultiCircle[nb][ils1+2].y;
        pp2.x = MultiCircle[nb][ils1].x;
        pp2.y = MultiCircle[nb][ils1].y;
        pp3.x = MultiCircle[nb][ils1+1].x;
        pp3.y = MultiCircle[nb][ils1+1].y;

        rd1=getRadius(&pp1,&pp2,&pp3);

        pp1.x = MultiCircle[nb][ils1+1].x;
        pp1.y = MultiCircle[nb][ils1+1].y;
        pp2.x = MultiCircle[nb][ils1].x;
        pp2.y = MultiCircle[nb][ils1].y;
        pp3.x = MultiCircle[nb][ils1-1].x;
        pp3.y = MultiCircle[nb][ils1-1].y;

        rd2=getRadius(&pp1,&pp2,&pp3);

        return (0.5*(rd1+rd2));
    }
}

double getRadius(const PXYZ *pt1, const PXYZ *pt2,const PXYZ *pt3, PXYZ &centr)
{
    //要是三点共线或者接近共线这个中心可找不到 TODO
    double rd,dls1,dls2,dls3,dls4;
    dls1=(pt1->x-pt2->x)*(pt1->x-pt2->x)+(pt1->y-pt2->y)*(pt1->y-pt2->y);
    dls2=(pt1->x-pt3->x)*(pt1->x-pt3->x)+(pt1->y-pt3->y)*(pt1->y-pt3->y);
    dls3=(pt2->x-pt3->x)*(pt2->x-pt3->x)+(pt2->y-pt3->y)*(pt2->y-pt3->y);
    dls4=2.0*(pt3->x * (pt2->y - pt1->y) + pt2->x * (pt1->y - pt3->y) 
        + pt1->x * (pt3->y-pt2->y));
    rd = sqrt(dls1*dls2*dls3)/(fabs(dls4)+1.0e-18);

    double xc12,yc12;
    xc12=0.5*(pt1->x+pt2->x);
    yc12=0.5*(pt1->y+pt2->y);

    double dis1c, dis;
    dis1c=sqrt((pt1->x-xc12)*(pt1->x-xc12)+(pt1->y-yc12)*(pt1->y-yc12));
    dis=sqrt(rd*rd-dis1c*dis1c);

    //   PXYZ tdir;
    PXYZ nml, cen1;       
    //  tdir.x=(pt1->x-xc12)/dis1c;
    // tdir.y=(pt1->y-yc12)/dis1c;       

    nml.x=(pt1->y-yc12)/dis1c;  
    nml.y=-(pt1->x-xc12)/dis1c;

    cen1.x= xc12+nml.x*dis;
    cen1.y= yc12+nml.y*dis;

    if(fabs(rd*rd-(cen1.x-pt3->x)*(cen1.x-pt3->x)-(cen1.y-pt3->y)*(cen1.y-pt3->y))<ERRS) {
        centr=cen1;
    }
    else {
        centr.x= xc12-nml.x*dis;
        centr.y= yc12-nml.y*dis;        
    }  

    double a1,b1,a2,b2;
    a1=pt2->x-pt1->x;
    b1=pt2->y-pt1->y;

    a2=pt3->x-pt1->x;
    b2=pt3->y-pt1->y;
    if(fabs(a1*b2-b1*a2)<1.0e-12) rd=1.0e15;
    return rd;
}
//double getRadius(const PXYZ *pt1, const PXYZ *pt2,const PXYZ *pt3, PXYZ &centr);	
double getWallCurvatureRadius(const PXYZ *pt, PXYZ &ndir)
{
    int ils1;
    double rd1=0.0,rd2=0.0, dis1,dis2;
    PXYZ pp1,pp2,pp3, cen1,cen2, ndir1,ndir2;
    int nb;
    ils1=getIndexAtWallPoint(pt,nb);
    if(ils1==NWallPts[nb]-1)
    {
        pp1.x = MultiCircle[nb][0].x;
        pp1.y = MultiCircle[nb][0].y;
        pp2.x = MultiCircle[nb][ils1].x;
        pp2.y = MultiCircle[nb][ils1].y;
        pp3.x = MultiCircle[nb][ils1-1].x;
        pp3.y = MultiCircle[nb][ils1-1].y;

        rd1=getRadius(&pp1,&pp2,&pp3,cen1);
        dis1=sqrt((pt->x-cen1.x)*(pt->x-cen1.x)+(pt->y-cen1.y)*(pt->y-cen1.y));
        ndir.x=(pt->x-cen1.x)/dis1;
        ndir.y=(pt->y-cen1.y)/dis1;		

        return rd1;
    }
    else if(ils1==NWallPts[nb]-2)
    {
        pp1.x = MultiCircle[nb][0].x;
        pp1.y = MultiCircle[nb][0].y;
        pp2.x = MultiCircle[nb][ils1].x;
        pp2.y = MultiCircle[nb][ils1].y;
        pp3.x = MultiCircle[nb][ils1+1].x;
        pp3.y = MultiCircle[nb][ils1+1].y;

        rd1=getRadius(&pp1,&pp2,&pp3,cen1);
        dis1=sqrt((pt->x-cen1.x)*(pt->x-cen1.x)+(pt->y-cen1.y)*(pt->y-cen1.y));
        ndir1.x=(pt->x-cen1.x)/dis1;
        ndir1.y=(pt->y-cen1.y)/dis1;	

        pp1.x = MultiCircle[nb][ils1+1].x;
        pp1.y = MultiCircle[nb][ils1+1].y;
        pp2.x = MultiCircle[nb][ils1].x;
        pp2.y = MultiCircle[nb][ils1].y;
        pp3.x = MultiCircle[nb][ils1-1].x;
        pp3.y = MultiCircle[nb][ils1-1].y;

        rd2=getRadius(&pp1,&pp2,&pp3, cen2);

        dis2=sqrt((pt->x-cen2.x)*(pt->x-cen2.x)+(pt->y-cen2.y)*(pt->y-cen2.y));
        ndir2.x=(pt->x-cen2.x)/dis2;
        ndir2.y=(pt->y-cen2.y)/dis2;

        ndir.x=0.5*(ndir1.x+ndir2.x);
        ndir.y=0.5*(ndir1.y+ndir2.y);			

        return (0.5*(rd1+rd2));
    }
    else if(ils1==0)
    {
        pp1.x = MultiCircle[nb][0].x;
        pp1.y = MultiCircle[nb][0].y;
        pp2.x = MultiCircle[nb][1].x;
        pp2.y = MultiCircle[nb][1].y;
        pp3.x = MultiCircle[nb][2].x;
        pp3.y = MultiCircle[nb][2].y;
        //	rd1=getRadius(&pp1,&pp2,&pp3);
        rd1=getRadius(&pp1,&pp2,&pp3,cen1);
        dis1=sqrt((pt->x-cen1.x)*(pt->x-cen1.x)+(pt->y-cen1.y)*(pt->y-cen1.y));
        ndir.x=(pt->x-cen1.x)/dis1;
        ndir.y=(pt->y-cen1.y)/dis1;	

        return rd1;
    }
    else
    {
        pp1.x = MultiCircle[nb][ils1+2].x;
        pp1.y = MultiCircle[nb][ils1+2].y;
        pp2.x = MultiCircle[nb][ils1].x;
        pp2.y = MultiCircle[nb][ils1].y;
        pp3.x = MultiCircle[nb][ils1+1].x;
        pp3.y = MultiCircle[nb][ils1+1].y;

        rd1=getRadius(&pp1,&pp2,&pp3,cen1);
        dis1=sqrt((pt->x-cen1.x)*(pt->x-cen1.x)+(pt->y-cen1.y)*(pt->y-cen1.y));
        ndir1.x=(pt->x-cen1.x)/dis1;
        ndir1.y=(pt->y-cen1.y)/dis1;

        pp1.x = MultiCircle[nb][ils1+1].x;
        pp1.y = MultiCircle[nb][ils1+1].y;
        pp2.x = MultiCircle[nb][ils1].x;
        pp2.y = MultiCircle[nb][ils1].y;
        pp3.x = MultiCircle[nb][ils1-1].x;
        pp3.y = MultiCircle[nb][ils1-1].y;

        rd2=getRadius(&pp1,&pp2,&pp3, cen2);

        dis2=sqrt((pt->x-cen2.x)*(pt->x-cen2.x)+(pt->y-cen2.y)*(pt->y-cen2.y));
        ndir2.x=(pt->x-cen2.x)/dis2;
        ndir2.y=(pt->y-cen2.y)/dis2;

        ndir.x=0.5*(ndir1.x+ndir2.x);
        ndir.y=0.5*(ndir1.y+ndir2.y);

        return (0.5*(rd1+rd2));
    }
}

//which leaf cell will the point locate in?
OctCell *findPointInCell(const PXYZ *point, OctCell *cellflow)
{
    int m;
    //double plevel,hc;
    double hcx, hcy;
    OctCell *incell = NULL;

    if(cellflow->reflag == 0)	
    {	          	    
        //	plevel=cellflow->level;		           
        //	hc=h/pow(2.0,plevel);
        hcx=cellflow->dx;
        hcy=cellflow->dy;

        if(point->x < cellflow->xc1 + 0.5 * hcx + ERRS
            && point->x > cellflow->xc1 - 0.5 * hcx - ERRS
            && point->y < cellflow->yc1 + 0.5 * hcy + ERRS
            && point->y > cellflow->yc1 - 0.5 * hcy - ERRS)				
        {					
            incell=cellflow;   
       //     return incell;
        } 
      //  else
        //    return NULL;			
    }
    else
    {
        for(m=0; m<4; ++m)
        {
            incell = findPointInCell(point, cellflow->children[m]);
            if(incell!=NULL)break;//2008.2.22加
        }
    }

    return incell;
}

extern OctCell *bodygrid;
void find_sol_refp(OctCell *pp[],OctCell *parent);
//确定在虚拟点gstc:parent相应的流体中的点p_sm属于哪个pcell 
void getGridIndex(OctCell *parent00, const PXYZ *p_sm, ptrOctCell &pcell)
{
   OctCell *pp[8]={NULL};
  find_sol_refp(pp,bodygrid+parent00->NOparent);
/*
    OctCell *pp[24]={NULL}, *parent; 
    parent=bodygrid+parent00->NOparent;
    pp[0]=EastNeighbor(parent);
    pp[1]=NorthNeighbor(pp[0]);
    pp[2]=NorthNeighbor(parent);
    pp[3]=WestNeighbor(pp[2]);
    pp[4]=WestNeighbor(parent);
    pp[5]=SouthNeighbor(pp[4]);
    pp[6]=SouthNeighbor(parent);
    pp[7]=EastNeighbor(pp[6]);

    pp[8]=NorthNeighbor(pp[1]);

    pp[9]=WestNeighbor(pp[8]);
    pp[10]=WestNeighbor(pp[9]);
    pp[11]=WestNeighbor(pp[10]);
    pp[12]=WestNeighbor(pp[3]);
    pp[13]=WestNeighbor(pp[4]);
    pp[14]=WestNeighbor(pp[5]);
    pp[15]=SouthNeighbor(pp[14]);
    pp[16]=EastNeighbor(pp[15]);
    pp[17]=EastNeighbor(pp[16]);
    pp[18]=EastNeighbor(pp[17]);
    pp[19]=EastNeighbor(pp[18]);
    pp[20]=NorthNeighbor(pp[19]);
    pp[21]=NorthNeighbor(pp[20]);
    pp[22]=NorthNeighbor(pp[21]);
    pp[23]=NorthNeighbor(pp[22]);
    */
    /*
    //================
    pp[24]= NorthNeighbor(pp[23]);
    pp[25]=WestNeighbor(pp[24]);
    pp[26]=WestNeighbor(pp[25]);
    pp[27]=WestNeighbor(pp[26]);  
    pp[28]=WestNeighbor(pp[27]);  
    pp[29]=WestNeighbor(pp[28]);    

    pp[30]=WestNeighbor(pp[11]);
    pp[31]=WestNeighbor(pp[12]);
    pp[32]=WestNeighbor(pp[13]);  
    pp[33]=WestNeighbor(pp[14]);  
    pp[34]=WestNeighbor(pp[15]);    

    pp[35]=SouthNeighbor(pp[34]);

    pp[36]=SouthNeighbor(pp[15]);
    pp[37]=SouthNeighbor(pp[16]);
    pp[38]=SouthNeighbor(pp[17]);
    pp[39]=SouthNeighbor(pp[18]);
    pp[40]=SouthNeighbor(pp[19]);
    pp[41]=EastNeighbor(pp[40]);
    pp[42]=NorthNeighbor(pp[41]);
    pp[43]=NorthNeighbor(pp[42]);  
    pp[44]=NorthNeighbor(pp[43]);  
    pp[45]=NorthNeighbor(pp[44]);  
    pp[46]=NorthNeighbor(pp[45]); 
    pp[47]=NorthNeighbor(pp[46]);  */

    OctCell *smInCell=NULL;//对称点所在的cell

    smInCell=findPointInCell(p_sm, bodygrid+parent00->NOparent);
    if(smInCell==NULL){ 
        for(int ii=0; ii<8; ++ii) {
            if(pp[ii]!=NULL) {
                smInCell=findPointInCell(p_sm, pp[ii]);
                if(smInCell!=NULL) break;
            }
        } 
    }
#ifdef DEBUG  
    if(smInCell==NULL){
        cout<<"010wrong find symitric point cell\n";
        cout<<parent00->xc1<<","<<parent00->yc1<<endl;
        cout<<p_sm->x<<","<<p_sm->y<<endl;
    }

#endif   
    pcell= smInCell;
//    if(smInCell==NULL)pcell=parent; //2012/3/36
}

void getGridIndex3(OctCell *parent00, const PXYZ *p_sm, ptrOctCell &pcell)
{
 //  OctCell *pp[8]={NULL};
 
  // find_sol_refp(pp,bodygrid+parent->NOparent);

    OctCell *pp[48]={NULL}, *parent; 
    parent=bodygrid+parent00->NOparent;
    pp[0]=EastNeighbor(parent);
    pp[1]=NorthNeighbor(pp[0]);
    pp[2]=NorthNeighbor(parent);
    pp[3]=WestNeighbor(pp[2]);
    pp[4]=WestNeighbor(parent);
    pp[5]=SouthNeighbor(pp[4]);
    pp[6]=SouthNeighbor(parent);
    pp[7]=EastNeighbor(pp[6]);

    pp[8]=NorthNeighbor(pp[1]);

    pp[9]=WestNeighbor(pp[8]);
    pp[10]=WestNeighbor(pp[9]);
    pp[11]=WestNeighbor(pp[10]);
    pp[12]=WestNeighbor(pp[3]);
    pp[13]=WestNeighbor(pp[4]);
    pp[14]=WestNeighbor(pp[5]);
    pp[15]=SouthNeighbor(pp[14]);
    pp[16]=EastNeighbor(pp[15]);
    pp[17]=EastNeighbor(pp[16]);
    pp[18]=EastNeighbor(pp[17]);
    pp[19]=EastNeighbor(pp[18]);
    pp[20]=NorthNeighbor(pp[19]);
    pp[21]=NorthNeighbor(pp[20]);
    pp[22]=NorthNeighbor(pp[21]);
    pp[23]=NorthNeighbor(pp[22]);
    //================
    pp[24]= NorthNeighbor(pp[23]);
    pp[25]=WestNeighbor(pp[24]);
    pp[26]=WestNeighbor(pp[25]);
    pp[27]=WestNeighbor(pp[26]);  
    pp[28]=WestNeighbor(pp[27]);  
    pp[29]=WestNeighbor(pp[28]);    

    pp[30]=WestNeighbor(pp[11]);
    pp[31]=WestNeighbor(pp[12]);
    pp[32]=WestNeighbor(pp[13]);  
    pp[33]=WestNeighbor(pp[14]);  
    pp[34]=WestNeighbor(pp[15]);    

    pp[35]=SouthNeighbor(pp[34]);

    pp[36]=SouthNeighbor(pp[15]);
    pp[37]=SouthNeighbor(pp[16]);
    pp[38]=SouthNeighbor(pp[17]);
    pp[39]=SouthNeighbor(pp[18]);
    pp[40]=SouthNeighbor(pp[19]);
    pp[41]=EastNeighbor(pp[40]);
    pp[42]=NorthNeighbor(pp[41]);
    pp[43]=NorthNeighbor(pp[42]);  
    pp[44]=NorthNeighbor(pp[43]);  
    pp[45]=NorthNeighbor(pp[44]);  
    pp[46]=NorthNeighbor(pp[45]); 

    OctCell *smInCell=NULL;//对称点所在的cell

    smInCell=findPointInCell(p_sm, bodygrid+parent00->NOparent);
    if(smInCell==NULL){ 
        for(int ii=0; ii<48; ++ii) {
            if(pp[ii]!=NULL) {
                smInCell=findPointInCell(p_sm, pp[ii]);
                if(smInCell!=NULL) break;
            }
        } 
    }
#ifdef DEBUG  
    if(smInCell==NULL){
        cout<<"020wrong find symitric point cell\n";
        cout<<parent->xc1<<","<<parent->yc1<<endl;
        cout<<p_sm->x<<","<<p_sm->y<<endl;
    }

#endif   
    pcell= smInCell;
 //   if(smInCell==NULL)pcell=parent; //2012/3/36
}
void getGridIndex2(OctCell *parent00, const PXYZ *p_sm, ptrOctCell &pcell)
{
 //  OctCell *pp[8]={NULL};
 
  // find_sol_refp(pp,bodygrid+parent->NOparent);

    OctCell *pp[24]={NULL}, *parent; 
    parent=bodygrid+parent00->NOparent;
    pp[0]=EastNeighbor(parent);
    pp[1]=NorthNeighbor(pp[0]);
    pp[2]=NorthNeighbor(parent);
    pp[3]=WestNeighbor(pp[2]);
    pp[4]=WestNeighbor(parent);
    pp[5]=SouthNeighbor(pp[4]);
    pp[6]=SouthNeighbor(parent);
    pp[7]=EastNeighbor(pp[6]);

    pp[8]=NorthNeighbor(pp[1]);

    pp[9]=WestNeighbor(pp[8]);
    pp[10]=WestNeighbor(pp[9]);
    pp[11]=WestNeighbor(pp[10]);
    pp[12]=WestNeighbor(pp[3]);
    pp[13]=WestNeighbor(pp[4]);
    pp[14]=WestNeighbor(pp[5]);
    pp[15]=SouthNeighbor(pp[14]);
    pp[16]=EastNeighbor(pp[15]);
    pp[17]=EastNeighbor(pp[16]);
    pp[18]=EastNeighbor(pp[17]);
    pp[19]=EastNeighbor(pp[18]);
    pp[20]=NorthNeighbor(pp[19]);
    pp[21]=NorthNeighbor(pp[20]);
    pp[22]=NorthNeighbor(pp[21]);
    pp[23]=NorthNeighbor(pp[22]);
    //================
    /*
    pp[24]= NorthNeighbor(pp[23]);
    pp[25]=WestNeighbor(pp[24]);
    pp[26]=WestNeighbor(pp[25]);
    pp[27]=WestNeighbor(pp[26]);  
    pp[28]=WestNeighbor(pp[27]);  
    pp[29]=WestNeighbor(pp[28]);    

    pp[30]=WestNeighbor(pp[11]);
    pp[31]=WestNeighbor(pp[12]);
    pp[32]=WestNeighbor(pp[13]);  
    pp[33]=WestNeighbor(pp[14]);  
    pp[34]=WestNeighbor(pp[15]);    

    pp[35]=SouthNeighbor(pp[34]);

    pp[36]=SouthNeighbor(pp[15]);
    pp[37]=SouthNeighbor(pp[16]);
    pp[38]=SouthNeighbor(pp[17]);
    pp[39]=SouthNeighbor(pp[18]);
    pp[40]=SouthNeighbor(pp[19]);
    pp[41]=EastNeighbor(pp[40]);
    pp[42]=NorthNeighbor(pp[41]);
    pp[43]=NorthNeighbor(pp[42]);  
    pp[44]=NorthNeighbor(pp[43]);  
    pp[45]=NorthNeighbor(pp[44]);  
    pp[46]=NorthNeighbor(pp[45]); 
    */

    OctCell *smInCell=NULL;//对称点所在的cell

    smInCell=findPointInCell(p_sm, bodygrid+parent00->NOparent);
    if(smInCell==NULL){ 
        for(int ii=0; ii<24; ++ii) {
            if(pp[ii]!=NULL) {
                smInCell=findPointInCell(p_sm, pp[ii]);
                if(smInCell!=NULL) break;
            }
        } 
    }
#ifdef DEBUG  
    if(smInCell==NULL){
        cout<<"020wrong find symitric point cell\n";
        cout<<parent->xc1<<","<<parent->yc1<<endl;
        cout<<p_sm->x<<","<<p_sm->y<<endl;
    }

#endif   
    pcell= smInCell;
 //   if(smInCell==NULL)pcell=parent; //2012/3/36
}
void insert_cell(OctCell *pp, set<OctCell*> &grid)
{
    if(pp->reflag==0){
        grid.insert(pp);
    }
    else {
        for(int m=0; m<4;++m) 
            insert_cell(pp->children[m], grid);
    }
}
void findIntpltCell(OctCell * incell, const PXYZ &gstp, const PXYZ &smtricp, vector< OctCell * > &grid) 
{
    PXYZ drct, drct1;  
    OctCell *pp[8]={NULL};
    set<OctCell*> grid_set;

    grid.clear();
    pp[0]=EastNeighbor(incell);
    pp[1]=NorthNeighbor(pp[0]);
    pp[2]=NorthNeighbor(incell);
    pp[3]=WestNeighbor(pp[2]);
    pp[4]=WestNeighbor(incell);
    pp[5]=SouthNeighbor(pp[4]);
    pp[6]=SouthNeighbor(incell);
    pp[7]=EastNeighbor(pp[6]);

    drct.x=smtricp.x-gstp.x;
    drct.y=smtricp.y-gstp.y;
    // cout<<iy<<" "<<ix<<endl;
    //cout<<smtricp.x<<" "<<smtricp.y<<endl;
    insert_cell(incell, grid_set);
    for(int ii=0;ii<8;++ii) insert_cell(pp[ii], grid_set);

    for(set<OctCell*>::iterator itset=grid_set.begin(); itset!=grid_set.end(); ++itset) {
        if((*itset)->flag==0) {
            drct1.x= (*itset)->xc1-gstp.x;
            drct1.y= (*itset)->yc1-gstp.y;  
            if(drct.x*drct1.x+drct.y*drct1.y+1.0 >= 1.0) {
                grid.push_back(*itset); 
            }          
        }
    }

#ifdef DEBUG
    int itm;
    if(grid.size()==0){
        cout << "wrong findIntpltCell"<< endl;
        cout<<gstp.x<<" "<<gstp.y<<", "<<smtricp.x<<"  "<<smtricp.y<<endl;
        cin>>itm;
    }
#endif
}   

//incell为反射点属于的控制体
/*
void findIntpltCell(OctCell * incell, const PXYZ &gstp, const PXYZ &smtricp, vector< OctCell * > &grid4c,int j) 
{
    PXYZ drct, drct1,pct;  
    OctCell *pp[8]={NULL};
    set<OctCell*> grid_set;
    //vector<OctCell *> grid;
    multimap<double, OctCell *> mt_map;
    double rtmp;
    int itmp=0;

    grid4c.clear();
    pp[0]=EastNeighbor(incell);
    pp[1]=NorthNeighbor(pp[0]);
    pp[2]=NorthNeighbor(incell);
    pp[3]=WestNeighbor(pp[2]);
    pp[4]=WestNeighbor(incell);
    pp[5]=SouthNeighbor(pp[4]);
    pp[6]=SouthNeighbor(incell);
    pp[7]=EastNeighbor(pp[6]);

    drct.x=smtricp.x-gstp.x;
    drct.y=smtricp.y-gstp.y;
    // cout<<iy<<" "<<ix<<endl;
    //cout<<smtricp.x<<" "<<smtricp.y<<endl;
    insert_cell(incell, grid_set);
    for(int ii=0;ii<8;++ii) insert_cell(pp[ii], grid_set);

    for(set<OctCell*>::iterator itset=grid_set.begin(); itset!=grid_set.end(); ++itset) {
        if((*itset)->flag==0) {
            pct.x=(*itset)->xc1;
            pct.y=(*itset)->yc1;
            drct1.x= pct.x-gstp.x;
            drct1.y= pct.y-gstp.y;  
            if(drct.x*drct1.x+drct.y*drct1.y+1.0 >= 1.0) {
               // grid.push_back(*itset); 
               rtmp=distance2p(smtricp,pct);
               mt_map.insert(make_pair(rtmp+40.0,*itset));
            }          
        }
    }

    for(multimap<double, OctCell*> ::iterator ite=mt_map.begin(); ite!=mt_map.end();++ite ) {
       // ntmp=mt_map.count(ite->first);
        //for(int cnt=0; cnt!=ntmp; ++cnt, ++ite)
        grid4c.push_back(ite->second);
        if(grid4c.size()==4)break;
    }
#ifdef DEBUG
    if(grid4c.size()!=4) {
        cout<<grid4c.size()<<"wrong find interpolation cell <4 ";
        cout<<grid4c[0]->xc1<<", "<<grid4c[0]->yc1<<"\n";
    }

    if(grid4c.size()==0) cout << "wrong findIntpltCell"<< endl;
#endif
}*/   

void gem(double *a, int M, double *u)
{	
    int i,j,k1,p,mm,flag=1;
    double eita,ui3;

    const double nrrs=1e-8;
    //double *r=(double *)malloc((M+1)*sizeof(double));
    double *r=new double[M+1];

    //-----M*(M+1)增广矩阵列主元Gauss消去--------/
    for(k1=0;k1<=M-2;k1++)            
    {
        ui3=0.0;

        for(mm=k1;mm<=M-1;mm++)
        {
            if(ui3<fabs(a[mm*(M+1)+k1]))
            {
                ui3=fabs(a[mm*(M+1)+k1]);
                p=mm;
            }
        }

        if(ui3 < nrrs) 
        {
            // printf("matrix singularity.\n");
            cout<<"matrix singularity.\n";
            flag=0;
            break;
        }
        if(p!=k1)
        {
            for(j=k1;j<=M;j++) 
            {
                ui3=a[k1*(M+1)+j];
                a[k1*(M+1)+j]=a[p*(M+1)+j];
                a[p*(M+1)+j]=ui3;
            }  
        }
        for(j=k1+1;j<=M;j++)
            r[j]=a[k1*(M+1)+j];
        for(i=k1+1;i<=M-1;i++)
        {
            eita=a[i*(M+1)+k1]/a[k1*(M+1)+k1];
            a[i*(M+1)+k1]=eita;
            for(j=k1+1;j<=M;j++)
                a[i*(M+1)+j]=a[i*(M+1)+j]-eita*r[j];
        }
    }
    //-----回代得到解向量-------
    if(flag)
    {
        for(i=M-1;i>=0;i--)
        {	
            u[i]=a[i*(M+1)+M];
            for(j=i+1;j<=M-1;j++)
                u[i]=u[i]-a[i*(M+1)+j]*u[j];
            u[i]=u[i]/a[i*(M+1)+i];
        }
    }
    //free(r);
    delete [] r;
}

void bilinearIntplt(vector< OctCell * > &grid, const PXYZ *px, JBBL &ingrid1)
{
    double amatrix[20];
    double ujie[4];
    SHBL shu;
    PXYZ pxdc; //px所对应的对称点坐标
    int i;

    for( i=0; i<4; i++)
    {
        amatrix[i * 5 + 0] =1.0;
        amatrix[i * 5 + 1] = grid[i]->xc1;
        amatrix[i * 5 + 2] = grid[i]->yc1;
        amatrix[i * 5 + 3] = grid[i]->xc1 * grid[i]->yc1;
    }

    for(i=0; i<4; i++)
    {
        amatrix[i * 5 + 4] = grid[i]->dof[0][0];
    }

    gem(amatrix, 4, ujie);

    shu.q = ujie[0] + ujie[1] * px->x + ujie[2] * px->y + ujie[3] * px->x * px->y;

    for( i=0; i<4; i++)
    {
        amatrix[i * 5 + 0] =1.0;
        amatrix[i * 5 + 1] = grid[i]->xc1;
        amatrix[i * 5 + 2] = grid[i]->yc1;
        amatrix[i * 5 + 3] = grid[i]->xc1 * grid[i]->yc1;
    }

    for(i=0; i<4; i++)
    {
        amatrix[i * 5 + 4] = grid[i]->dof[1][0];
    }

    gem(amatrix, 4, ujie);

    shu.qu = ujie[0] + ujie[1] * px->x + ujie[2] * px->y + ujie[3] * px->x * px->y;

    for( i=0; i<4; i++)
    {
        amatrix[i * 5 + 0] =1.0;
        amatrix[i * 5 + 1] = grid[i]->xc1;
        amatrix[i * 5 + 2] = grid[i]->yc1;
        amatrix[i * 5 + 3] = grid[i]->xc1 * grid[i]->yc1;
    }

    for(i=0; i<4; i++)
    {
        amatrix[i * 5 + 4] = grid[i]->dof[2][0];
    }

    gem(amatrix, 4, ujie);

    shu.qv = ujie[0] + ujie[1] * px->x + ujie[2] * px->y + ujie[3] * px->x * px->y;

    for( i=0; i<4; i++)
    {
        amatrix[i * 5 + 0] =1.0;
        amatrix[i * 5 + 1] = grid[i]->xc1;
        amatrix[i * 5 + 2] = grid[i]->yc1;
        amatrix[i * 5 + 3] = grid[i]->xc1 * grid[i]->yc1;
    }

    for(i=0; i<4; i++)
    {
        amatrix[i * 5 + 4] = grid[i]->dof[3][0];
    }

    gem(amatrix, 4, ujie);

    shu.te = ujie[0] + ujie[1] * px->x + ujie[2] * px->y + ujie[3] * px->x * px->y;

    shblToJbbl(&shu, &ingrid1);
}

//得到基本变量结果  //是否换成对基本变量进行插值比较好 TODO
void chazhi4smtrictPshbl(vector< OctCell * > &grid, const PXYZ &smtrict, JBBL &jbu )
{
    double distance;
    vector<double> expMinusR;
    double sumR = 0.0; 
    int nmc1=grid.size();  

    SHBL shu;
    for(int i=0; i < nmc1; ++i) {
        distance=sqrt( (smtrict.x-grid[i]->xc1) * (smtrict.x-grid[i]->xc1) 
            + (smtrict.y-grid[i]->yc1) * (smtrict.y-grid[i]->yc1));
        // expMinusR.push_back(exp(-distance));
        expMinusR.push_back(1.0/(1.0e-18+distance));
        sumR += expMinusR[i];

    }

    for(int i = 0; i < nmc1; ++i)
    {
        shu.q += expMinusR[i] * grid[i]->dof[0][0]/sumR;
        shu.qu += expMinusR[i] * grid[i]->dof[1][0]/sumR;
        shu.qv += expMinusR[i] * grid[i]->dof[2][0]/sumR;
        shu.te += expMinusR[i] * grid[i]->dof[3][0]/sumR;
    }

    shblToJbbl(&shu,&jbu);   
}

void getGhostValueMagazine(OctCell * ghostcell, int lyx)
{
    PXYZ nml1;
    PXYZ symtrcpt,jiaodian,wnph; 
    PXYZ ghostpt, plsA, plsB;

    int jj,kk,m,numbergrid,lsflag,lsint, ils;
    int nb;
    //Rieman problem界面量
    //    double dil,dir,ui,prei;
    JBBL ingrid1,jbu1; 

    JBBL jbgst1;

    SHBL shgst1;

    //	OctCell *leafp;
    OctCell *pincell=NULL;
    vector< OctCell * > grid;

    double pwh1, vnsolid,pw1;	
    double  hc;

    ghostpt.x = ghostcell->xc1;
    ghostpt.y = ghostcell->yc1;	

    //hc=max2(ghostcell->dx, ghostcell->dy);
    hc=sqrt(ghostcell->dx*ghostcell->dx+ghostcell->dy*ghostcell->dy);

    ils = getIndexAtWallPoint(&ghostpt,nb); //TODO

    if(ils == NWallPts[nb]-1) 
    {
        plsA.x = MultiCircle[nb][ils].x; 
        plsA.y = MultiCircle[nb][ils].y;
        plsB.x = MultiCircle[nb][0].x; 
        plsB.y = MultiCircle[nb][0].y;
    }
    else
    {
        plsA.x = MultiCircle[nb][ils].x; 
        plsA.y = MultiCircle[nb][ils].y;
        plsB.x = MultiCircle[nb][ils+1].x; 
        plsB.y = MultiCircle[nb][ils+1].y;
    }

    getsmtrctAB(&ghostpt, &plsA, &plsB,
        &symtrcpt, &jiaodian, &nml1);

#ifdef DEBUG
    //cout<<"("<<ghostpt.x<<", "<<ghostpt.y<<")  ("<<jiaodian.x<<", "<<jiaodian.y<<")   ("
    // << symtrcpt.x<<", "<<symtrcpt.y<<")\n";
#endif

    //void getGridIndex(OctCell *gstc, const PXYZ *p_sm, int &iy, int &ix)
    //    getGridIndex(ghostcell, &symtrcpt, jj, kk);
    //void getGridIndex(OctCell *parent, const PXYZ *p_sm, ptrOctCell &pcell)

    getGridIndex(ghostcell, &symtrcpt, pincell);

    //----------------------------------------------------------begin
    //void findIntpltCell(int iy, int ix, const PXYZ &gstp, const PXYZ &smtricp, vector< OctCell * > &grid) 

    //   findIntpltCell(jj, kk, ghostpt, symtrcpt, grid);
    //incell为反射点属于的控制体
    //void findIntpltCell(const OctCell * incell, const PXYZ &gstp, const PXYZ &smtricp, vector< OctCell * > &grid)     
    findIntpltCell(pincell, ghostpt, symtrcpt, grid);

    // 	numbergrid = grid.size();

    //cout<<"numbergrid"<<numbergrid<<endl;

    // 	void chazhi4smtrictPshbl(vector< OctCell * > &grid, const PXYZ &smtrict, JBBL &jbu )

    chazhi4smtrictPshbl(grid, symtrcpt, ingrid1); 

    //2011-9.9	  									
    //    	InterpolateForMoving(grid, &symtrcpt, numbergrid, &ingrid1); 
    //------------------------------------------------------end-----------

    //	linearInterpolate(grid, &symtrcpt, &ingrid1);  

    reflection(&ingrid1, &nml1, &jbgst1);
    //////////////////				
    jbblToShbl(&jbgst1, &shgst1);

    ghostcell->dof[0][0] = shgst1.q;
    ghostcell->dof[1][0] = shgst1.qu;
    ghostcell->dof[2][0]  = shgst1.qv;
    ghostcell->dof[3][0] = shgst1.te;
    for(int itm=0; itm<4; ++itm) {
        for(int jtm=1;jtm<nDOF; ++jtm) {
            ghostcell->dof[itm][jtm]=0.0;
        }
    }
}



void getGhostValueMagazine(OctCell * ghostcell)
{
    PXYZ nml1;
    PXYZ symtrcpt,jiaodian,wnph; 
    PXYZ ghostpt, plsA, plsB;

    int jj,kk,m,numbergrid,lsflag,lsint, ils;
    int nb;
//Rieman problem界面量
//    double dil,dir,ui,prei;
    JBBL ingrid1,jbu1; 

    JBBL jbgst1;

    SHBL shgst1;

    //	OctCell *leafp;
    OctCell *pincell=NULL;
    vector< OctCell * > grid;

    double pwh1, vnsolid,pw1;	
    double  hc;

    ghostpt.x = ghostcell->xc1;
    ghostpt.y = ghostcell->yc1;	

    //hc=max2(ghostcell->dx, ghostcell->dy);
    hc=sqrt(ghostcell->dx*ghostcell->dx+ghostcell->dy*ghostcell->dy);

    if(qvlvModify == 0)
    {
        ils = getIndexAtWallPoint(&ghostpt,nb); //TODO

        if(ils == NWallPts[nb]-1) 
        {
            plsA.x = MultiCircle[nb][ils].x; 
            plsA.y = MultiCircle[nb][ils].y;
            plsB.x = MultiCircle[nb][0].x; 
            plsB.y = MultiCircle[nb][0].y;
        }
        else
        {
            plsA.x = MultiCircle[nb][ils].x; 
            plsA.y = MultiCircle[nb][ils].y;
            plsB.x = MultiCircle[nb][ils+1].x; 
            plsB.y = MultiCircle[nb][ils+1].y;
        }

        getsmtrctAB(&ghostpt, &plsA, &plsB,
            &symtrcpt, &jiaodian, &nml1);

#ifdef DEBUG
        //cout<<"("<<ghostpt.x<<", "<<ghostpt.y<<")  ("<<jiaodian.x<<", "<<jiaodian.y<<")   ("
        // << symtrcpt.x<<", "<<symtrcpt.y<<")\n";
#endif

        //void getGridIndex(OctCell *gstc, const PXYZ *p_sm, int &iy, int &ix)
        //    getGridIndex(ghostcell, &symtrcpt, jj, kk);
        //void getGridIndex(OctCell *parent, const PXYZ *p_sm, ptrOctCell &pcell)

        getGridIndex(ghostcell, &symtrcpt, pincell);

        //----------------------------------------------------------begin
        //void findIntpltCell(int iy, int ix, const PXYZ &gstp, const PXYZ &smtricp, vector< OctCell * > &grid) 

        //   findIntpltCell(jj, kk, ghostpt, symtrcpt, grid);
        //incell为反射点属于的控制体
        //void findIntpltCell(const OctCell * incell, const PXYZ &gstp, const PXYZ &smtricp, vector< OctCell * > &grid)     
        findIntpltCell(pincell, ghostpt, symtrcpt, grid);

        // 	numbergrid = grid.size();

        //cout<<"numbergrid"<<numbergrid<<endl;

        // 	void chazhi4smtrictPshbl(vector< OctCell * > &grid, const PXYZ &smtrict, JBBL &jbu )

        chazhi4smtrictPshbl(grid, symtrcpt, ingrid1); 

        //2011-9.9	  									
        //    	InterpolateForMoving(grid, &symtrcpt, numbergrid, &ingrid1); 
        //------------------------------------------------------end-----------

        //	linearInterpolate(grid, &symtrcpt, &ingrid1);  

        reflection(&ingrid1, &nml1, &jbgst1);
    }//------------------------------------------------------------------------
    else if(qvlvModify == 1)
    {
        ils = getIndexAtWallPoint(&ghostpt,nb); //TODO

        if(ils == NWallPts[nb]-1) 
        {
            plsA.x = MultiCircle[nb][ils].x; 
            plsA.y = MultiCircle[nb][ils].y;
            plsB.x = MultiCircle[nb][0].x; 
            plsB.y = MultiCircle[nb][0].y;
        }
        else
        {
            plsA.x = MultiCircle[nb][ils].x; 
            plsA.y = MultiCircle[nb][ils].y;
            plsB.x = MultiCircle[nb][ils+1].x; 
            plsB.y = MultiCircle[nb][ils+1].y;
        }
        /////////////////////////////////////2008-12-6//////////////////////////////////////////////

        getsmtrctAB(&ghostpt, &plsA, &plsB,
            &symtrcpt, &jiaodian, &nml1);
        //确定流体中的点p_sm属于哪个cell

        //void getGridIndex(OctCell *gstc, const PXYZ *p_sm, int &iy, int &ix)
        //	getGridIndex(ghostcell, &symtrcpt, jj, kk);
        getGridIndex(ghostcell, &symtrcpt, pincell);       	

        //  findIntpltCell(jj, kk, ghostpt, symtrcpt, grid);
        findIntpltCell(pincell, ghostpt, symtrcpt, grid);
        //----------------------------------------------------------begin
        //	numbergrid = grid.size();
        // 	void chazhi4smtrictPshbl(vector< OctCell * > &grid, const PXYZ &smtrict, JBBL &jbu )

        chazhi4smtrictPshbl(grid, symtrcpt, ingrid1); 

        // 	InterpolateForMoving(grid, &symtrcpt, numbergrid, &ingrid1); 
        //------------------------------------------------------end--------
        vnsolid=0.0;						
        reflectionForMoving6(hc,&ingrid1, &jiaodian, &ghostpt, 
            &nml1, vnsolid, &jbgst1);
//void  solveARS(double dl,double dr,double gl,double gr,double bl,
//			   double br,double ul,double ur,double pl,double pr,
//			   double *dil,double *dir,double *ui,double *pi);
        //////////////////////////2008-12-6//////2008-12-6///////////////////////////////////////////
    }
    //////////////////				
    jbblToShbl(&jbgst1, &shgst1);

    ghostcell->dof[0][0] = shgst1.q;
    ghostcell->dof[1][0] = shgst1.qu;
    ghostcell->dof[2][0]  = shgst1.qv;
    ghostcell->dof[3][0] = shgst1.te;
    for(int itm=0; itm<4; ++itm) {
        for(int jtm=1;jtm<nDOF; ++jtm) {
            ghostcell->dof[itm][jtm]=0.0;
        }
    }
}



void wall_treatment()
{

    BndryNode *current;
    current = HeadFaceCell;
    while(current != NULL)
    {
        getGhostValueMagazine(current->cell);
        current->cell->bAvg=true;  //2012/4/27 add
        current = current->next;
    }
}

//extern Vec2D ExtPnt0;
//extern Vec2D ExtPnt1;
//extern Vec2D ExtPnt2;
//extern Vec2D ExtPnt3;

void getGhostValueExtWall(OctCell * ghostcell)
{
    PXYZ nml1;
    PXYZ symtrcpt,jiaodian,wnph; 
    PXYZ ghostpt, plsA, plsB;

    int  ils;
    //Rieman problem界面量
    //    double dil,dir,ui,prei;
    JBBL ingrid1,jbu1; 

    JBBL jbgst1;

    SHBL shgst1;

    //	OctCell *leafp;
    OctCell *pincell=NULL;
    vector< OctCell * > grid;

//    double pwh1, vnsolid,pw1;	
 //   double  hc;

    ghostpt.x = ghostcell->xc1;
    ghostpt.y = ghostcell->yc1;	

    //hc=max2(ghostcell->dx, ghostcell->dy);
//    hc=sqrt(ghostcell->dx*ghostcell->dx+ghostcell->dy*ghostcell->dy);

    ils = getIndexAtWallPoint(&ghostpt, NPT_EXT_WALL, ext_wall); //TODO

    if(ils == NPT_EXT_WALL-1) 
    {
        plsA.x = ext_wall[ils].x; 
        plsA.y = ext_wall[ils].y;
        plsB.x = ext_wall[0].x; 
        plsB.y = ext_wall[0].y;
    }
    else
    {
        plsA.x = ext_wall[ils].x; 
        plsA.y = ext_wall[ils].y;
        plsB.x = ext_wall[ils+1].x; 
        plsB.y = ext_wall[ils+1].y;
    }

    getsmtrctAB(&ghostpt, &plsA, &plsB,
        &symtrcpt, &jiaodian, &nml1);

#ifdef DEBUG
    //cout<<"("<<ghostpt.x<<", "<<ghostpt.y<<")  ("<<jiaodian.x<<", "<<jiaodian.y<<")   ("
    // << symtrcpt.x<<", "<<symtrcpt.y<<")\n";
#endif

    //void getGridIndex(OctCell *gstc, const PXYZ *p_sm, int &iy, int &ix)
    //    getGridIndex(ghostcell, &symtrcpt, jj, kk);
    //void getGridIndex(OctCell *parent, const PXYZ *p_sm, ptrOctCell &pcell)

    getGridIndex2(ghostcell, &symtrcpt, pincell);

    //----------------------------------------------------------begin
    //void findIntpltCell(int iy, int ix, const PXYZ &gstp, const PXYZ &smtricp, vector< OctCell * > &grid) 

    //   findIntpltCell(jj, kk, ghostpt, symtrcpt, grid);
    //incell为反射点属于的控制体
    //void findIntpltCell(const OctCell * incell, const PXYZ &gstp, const PXYZ &smtricp, vector< OctCell * > &grid)     
    findIntpltCell(pincell, ghostpt, symtrcpt, grid);

    // 	numbergrid = grid.size();

    //cout<<"numbergrid"<<numbergrid<<endl;

    // 	void chazhi4smtrictPshbl(vector< OctCell * > &grid, const PXYZ &smtrict, JBBL &jbu )

    chazhi4smtrictPshbl(grid, symtrcpt, ingrid1); 

    //2011-9.9	  									
    //    	InterpolateForMoving(grid, &symtrcpt, numbergrid, &ingrid1); 
    //------------------------------------------------------end-----------

    //	linearInterpolate(grid, &symtrcpt, &ingrid1);  

    reflection(&ingrid1, &nml1, &jbgst1);
    //////////////////				
    jbblToShbl(&jbgst1, &shgst1);

    ghostcell->dof[0][0] = shgst1.q;
    ghostcell->dof[1][0] = shgst1.qu;
    ghostcell->dof[2][0]  = shgst1.qv;
    ghostcell->dof[3][0] = shgst1.te;
    for(int itm=0; itm<4; ++itm) {
        for(int jtm=1;jtm<nDOF; ++jtm) {
            ghostcell->dof[itm][jtm]=0.0;
        }
    }
}

//extern vector<OctCell*> ExtWallBndry;
void setExtWallBndry()
{
    for(vector<OctCell*>::size_type isz=0; isz<ExtWallBndry.size(); ++isz) {
        getGhostValueExtWall(ExtWallBndry[isz]);
    }

}
//int newFlagpoint(double x, double y);

double distance2p(const PXYZ &p1, const PXYZ &p2);

//void rbf_mq( int np, PXYZ *pxy, PXYZ &xy, double *u, double c, double &u_sol);
//void rbf_mqUxUy( int np, PXYZ *pxy, PXYZ &xy, double *u, double c, 
  //  double &u_sol, double &Ux, double &Uy);
//extern vector<OctCell*> inter_cell;
//intersection cell, wall point in this cell
static set<OctCell*> interface_cell_set; //可能的,不一定就是

void insert_cell2(OctCell *pp, set<OctCell*> &grid);
//extern BndryNode *HeadFaceCell;
void form_inter_cell_set()
{
    BndryNode *current;
    OctCell *lsbl;
    current =HeadFaceCell;
    interface_cell_set.clear();
    OctCell *pp[8]={NULL};
    while(current != NULL)
    {
        lsbl=current->cell;
        find_sol_refp(pp,bodygrid+lsbl->NOparent);
        insert_cell2(bodygrid+lsbl->NOparent,interface_cell_set);
        for(int i=0;i<8;++i)
            insert_cell2(pp[i],interface_cell_set);

        current = current->next;
    }
}

void wherePoint(PXYZ *xy, ptrOctCell &incell)
{
    double xc1,yc1, dx1,dy1;

    OctCell *lsbl;
    for(set<OctCell*>::iterator itset=interface_cell_set.begin();
        itset!=interface_cell_set.end();++itset){
        lsbl =*itset;

        xc1=lsbl->xc1;
        yc1=lsbl->yc1;
        dx1=lsbl->dx;
        dy1=lsbl->dy;

        if(xy->x >= xc1-0.5*dx1-1e-10 && xy->x <= xc1+0.5*dx1+1e-10 
            && xy->y >= yc1-0.5*dy1-1e-10 && xy->y <= yc1+0.5*dy1+1e-10)
        {
            incell=lsbl;
            break;
        }
    }
}

void wherePointAll()
{
    int i;

    for (int jj=0; jj<N_BODY;++jj){
        for(i=0; i < NWallPts[jj]; i++)
        {
            wherePoint(&MultiCircle[jj][i],MultiCompany[jj][i].incell );
        }
        MultiCompany[jj][NWallPts[jj]]=MultiCompany[jj][0];
    }
}

void computeNormalArea()
{
    int i=0;

    for (int jj=0; jj<N_BODY;++jj){
        for(i=0; i<NWallPts[jj]; i++)
        {
            MultiCompany[jj][i].normS[0]= MultiCircle[jj][i+1].y-MultiCircle[jj][i].y;
            MultiCompany[jj][i].normS[1]= MultiCircle[jj][i].x-MultiCircle[jj][i+1].x;
        }
    }
}

//extern double CL,Cd;
extern double sumVorM;
extern string str0;
void computeCL_Cd(int ni)
{
    int i;
    ofstream clcd_f;
    string str_T;
    str_T=str0+"CL_CD.out";
    if(ni==0) clcd_f.open(str_T.c_str());
    else clcd_f.open(str_T.c_str(), ofstream::app);
    /*
       double sum_fx=0.0;
       double sum_fy=0.0;

       for(i=0;i<NPOINTCIRCLE;i++)
       {
       sum_fx += Company[i].normS[0] * Company[i].pre;
       sum_fy += Company[i].normS[1] * Company[i].pre;
       }

       CL = 2.0*sum_fy;
       Cd = 2.0*sum_fx;
       */

    double fn[N_BODY+1]={0.0};
    double ft[N_BODY+1]={0.0};
    //	double angleA=AOA*PI/180.0;
    for(int wy=0;wy<N_BODY;++wy){
        for(i=0;i<NWallPts[wy];i++)
        {
            ft[wy] -= MultiCompany[wy][i].normS[0]
                * MultiCompany[wy][i].cp;
            fn[wy] -= MultiCompany[wy][i].normS[1] * MultiCompany[wy][i].cp;
        }
        //	CL =(fn*cos(angleA)-ft*sin(angleA));
        //  Cd =(fn*sin(angleA)+ft*cos(angleA));
        //	CL=fn[wy];
        //	Cd=ft[wy];
        cout<<"CL= "<<fn[wy]<<", Cd= "
            <<ft[wy]<<"  "<<"vor_sum= "<<sumVorM;
        clcd_f<<setw(15)<<fn[wy]<<setw(15)<<ft[wy]<<setw(15)<<sumVorM;
    }
    cout<<"\n";
    clcd_f<<"\n";
    clcd_f.close();
}



//void chazhi4smtrictPshbl(vector< OctCell * > &grid, const PXYZ &smtrict, JBBL &jbu )
void find_sol_4nbr(OctCell *pp[],OctCell *parent);
void computWallPresCp()
{
    OctCell *lsbl;
    OctCell *nbrcell[4]={NULL};
    //   OctCell *nbr[4]={NULL};
    vector<OctCell*> nbr;
    Vec2D dr;
    JBBL jbu;

    for (int jj=0; jj<N_BODY;++jj){
        for(int i=0; i<NWallPts[jj]; ++i){
            lsbl=MultiCompany[jj][i].incell;
            if(lsbl==NULL) {
                cout<<i<<" "<<MultiCircle[jj][i].x<<" "<<MultiCircle[jj][i].y<<"--";
                cout<<"wrong-computWallPresCp\n";break;}
            find_sol_4nbr(nbrcell,lsbl);
            nbr.clear();
            for(int inb=0;inb<4;++inb){
                if(nbrcell[inb]!=NULL ){
                    dr[0]=nbrcell[inb]->xc1-MultiCircle[jj][i].x;
                    dr[1]=nbrcell[inb]->yc1-MultiCircle[jj][i].y;
                    if(dr.dot(MultiCompany[jj][i].normS)>=0.0
                        && nbrcell[inb]->flag==0){
                        nbr.push_back(nbrcell[inb]);
                    }
                }
            }
            chazhi4smtrictPshbl(nbr,MultiCircle[jj][i],jbu);
            MultiCompany[jj][i].pre=jbu.p;
            MultiCompany[jj][i].cp=(jbu.p-initP)/(0.5*initQ*FreeMa*FreeMa);
        }
    }

}



//===============================================
/*
//=========================================
//r后面没有用
void findIntpltCell_cmplx(OctCell *cellofpoint, const PXYZ *point, vector<OctCell *>&grid);
//参考点不是镜像点
void getGhostValueMagazine_ref(OctCell * ghostcell)
{
    PXYZ nml1;
    PXYZ symtrcpt,jiaodian,wnph; 
    PXYZ ghostpt, plsA, plsB;

    int jj,kk,m,numbergrid,lsflag,lsint, ils;

    JBBL ingrid1,jbu1; 

    JBBL jbgst1;

    SHBL shgst1;

    //	OctCell *leafp;
    OctCell *pincell=NULL;
    vector< OctCell * > grid;

    double pwh1, vnsolid,pw1;	
    double  hr,hg;

    ghostpt.x = ghostcell->xc1;
    ghostpt.y = ghostcell->yc1;	

    hr=0.5*sqrt(ghostcell->dx*ghostcell->dx+ghostcell->dy*ghostcell->dy);

    if(qvlvModify == 0)
    {
        ils = getIndexAtWallPoint(&ghostpt); //TODO

        if(ils == NumberAirfoilPoint-1) 
        {
            plsA.x = wallpoint[ils].x; 
            plsA.y = wallpoint[ils].y;
            plsB.x = wallpoint[0].x; 
            plsB.y = wallpoint[0].y;
        }
        else
        {
            plsA.x = wallpoint[ils].x; 
            plsA.y = wallpoint[ils].y;
            plsB.x = wallpoint[ils+1].x; 
            plsB.y = wallpoint[ils+1].y;
        }

        getsmtrctAB(hr,&ghostpt, &plsA, &plsB,
            &symtrcpt, &jiaodian, &nml1,&hg);

#ifdef DEBUG
        //cout<<"("<<ghostpt.x<<", "<<ghostpt.y<<")  ("<<jiaodian.x<<", "<<jiaodian.y<<")   ("
        // << symtrcpt.x<<", "<<symtrcpt.y<<")\n";
#endif

        //void getGridIndex(OctCell *gstc, const PXYZ *p_sm, int &iy, int &ix)
        //    getGridIndex(ghostcell, &symtrcpt, jj, kk);
        //void getGridIndex(OctCell *parent, const PXYZ *p_sm, ptrOctCell &pcell)

        getGridIndex(ghostcell, &symtrcpt, pincell);

        //----------------------------------------------------------begin
        //void findIntpltCell(int iy, int ix, const PXYZ &gstp, const PXYZ &smtricp, vector< OctCell * > &grid) 

        //   findIntpltCell(jj, kk, ghostpt, symtrcpt, grid);
        //incell为反射点属于的控制体
        //void findIntpltCell(const OctCell * incell, const PXYZ &gstp, const PXYZ &smtricp, vector< OctCell * > &grid)     
        findIntpltCell(pincell, ghostpt, symtrcpt, grid,0);
       // findIntpltCell_cmplx(pincell, &symtrcpt, grid);

        // 	numbergrid = grid.size();

        //cout<<"numbergrid"<<numbergrid<<endl;

        // 	void chazhi4smtrictPshbl(vector< OctCell * > &grid, const PXYZ &smtrict, JBBL &jbu )

        chazhi4smtrictPshbl(grid, symtrcpt, ingrid1); 
       // cout<<ingrid1.q<<"?";
        //bilinearIntplt(grid, &symtrcpt, ingrid1); 
        //cout<<ingrid1.q<<" ";
       // for(int tt=0;tt<4;++tt)
         //   cout<<grid[tt]->xc1<<" "<<grid[tt]->yc1<<", ";
       // cout<<endl;
        //2011-9.9	  									
        //    	InterpolateForMoving(grid, &symtrcpt, numbergrid, &ingrid1); 
        //------------------------------------------------------end-----------

        //	linearInterpolate(grid, &symtrcpt, &ingrid1);  

        reflection(hr,hg,&ingrid1, &nml1, &jbgst1);
    }//------------------------------------------------------------------------
    else if(qvlvModify == 1)
    {
        /////////////////////////////////////2008-12-6//////////////////////////////////////////////
        ils = getIndexAtWallPoint(&ghostpt);

        if(ils == NumberAirfoilPoint-1) 
        {
            plsA.x = wallpoint[ils].x; 
            plsA.y = wallpoint[ils].y;
            plsB.x = wallpoint[0].x; 
            plsB.y = wallpoint[0].y;
        }
        else
        {
            plsA.x = wallpoint[ils].x; 
            plsA.y = wallpoint[ils].y;
            plsB.x = wallpoint[ils+1].x; 
            plsB.y = wallpoint[ils+1].y;
        }

        getsmtrctAB(hr,&ghostpt, &plsA, &plsB,
            &symtrcpt, &jiaodian, &nml1,&hg);
        //确定流体中的点p_sm属于哪个cell

        //void getGridIndex(OctCell *gstc, const PXYZ *p_sm, int &iy, int &ix)
        //	getGridIndex(ghostcell, &symtrcpt, jj, kk);
        getGridIndex(ghostcell, &symtrcpt, pincell);       	

        //  findIntpltCell(jj, kk, ghostpt, symtrcpt, grid);
        findIntpltCell(pincell, ghostpt, symtrcpt, grid,0);
        //findIntpltCell_cmplx(pincell, &symtrcpt, grid);
        //----------------------------------------------------------begin
        //	numbergrid = grid.size();
        // 	void chazhi4smtrictPshbl(vector< OctCell * > &grid, const PXYZ &smtrict, JBBL &jbu )

        chazhi4smtrictPshbl(grid, symtrcpt, ingrid1); 
       //bilinearIntplt(grid, &symtrcpt, ingrid1); 

        // 	InterpolateForMoving(grid, &symtrcpt, numbergrid, &ingrid1); 
        //------------------------------------------------------end--------
        vnsolid=0.0;						
        reflectionForMoving6(hr,hg,&ingrid1, &jiaodian, &ghostpt, 
            &nml1, vnsolid, &jbgst1);

        //////////////////////////2008-12-6//////2008-12-6///////////////////////////////////////////
    }
    //////////////////				
    jbblToShbl(&jbgst1, &shgst1);

    ghostcell->dof[0][0] = shgst1.q;
    ghostcell->dof[1][0] = shgst1.qu;
    ghostcell->dof[2][0]  = shgst1.qv;
    ghostcell->dof[3][0] = shgst1.te;
    for(int itm=0; itm<4; ++itm) {
        for(int jtm=1;jtm<nDOF; ++jtm) {
            ghostcell->dof[itm][jtm]=0.0;
        }
    }
}
*/
//++++++++++++++++++++++++++++++++++++后面没有了+++++++++
/*
void getCp(void)
{
    double Cpressure[NumberAirfoilPoint], press1;
    FILE *fpxx;

    int i, lsint, numbergrid, jj, kk, lsflag;
    static int ljm=0;

    PXYZ symtrcpt;
    OctCell *grid[4], *leafp;
    JBBL ingrid1;

    for(i = 0; i < NumberAirfoilPoint;i++)
    {
        symtrcpt.x = wallpoint[i].x;
        symtrcpt.y = wallpoint[i].y;
        getGridIndex(&symtrcpt, &POINTup0, &jj, &kk);

        leafp = findPointInCell(&symtrcpt, &bodygrid[(jj-1) * Nx + kk]);
        findIntplt8point(leafp, &symtrcpt, grid);
        //----------------------------------------------------------begin
        numbergrid = 0;
        if(grid[1] == NULL)
        {
            grid[numbergrid++] = grid[0];
        }
        else
        {
            for(lsint = 0; lsint < 4; lsint++)						
            {											
                lsflag=grid[lsint]->flag;
                if(lsflag == 0 || lsflag == 2 || lsflag == 4) grid[numbergrid++] = grid[lsint];					
            }			
        }			
        //	if(numbergrid<1)printf("dfdfdfd\n"), getch();
        InterpolateForMoving(grid, &symtrcpt, numbergrid, &ingrid1); 
        press1 = ingrid1.p;
        Cpressure[i]=(press1-initP)/(0.5*initQ*INITMach*INITMach);
    }

    if((fpxx = fopen(filename2[ljm],"w+")) == NULL)
    {
        printf("Can not open \n");
        exit(1);
    }

    ljm++;

    fprintf(fpxx,"TITLE = \"Cp UPPER\"\n");
    fprintf(fpxx,"VARIABLES = \"X\", \"Cp\"\n");
    fprintf(fpxx,"ZONE T= \"UPPER SURF.\", I = %d, F=POINT\n",
        NumberAirfoilPoint-HalfNumberAirfoilPoint+1);
    for(i = HalfNumberAirfoilPoint; i < NumberAirfoilPoint; i++)
        fprintf(fpxx,"%.8f  %.8f\n", wallpoint[i].x, Cpressure[i]);
    fprintf(fpxx,"%.8f  %.8f\n", wallpoint[0].x, Cpressure[0]);

    fprintf(fpxx,"TITLE = \"Cp LOWER\"\n");
    fprintf(fpxx,"VARIABLES = \"X\", \"Cp\"\n");
    fprintf(fpxx,"ZONE T= \"LOWER SURF.\", I = %d, F=POINT\n", HalfNumberAirfoilPoint+1);
    for(i = 0; i <= HalfNumberAirfoilPoint; i++)
        fprintf(fpxx,"%.8f  %.8f\n", wallpoint[i].x, Cpressure[i]);

    fclose(fpxx);
}
*/
//void computWallPresCp()
//{
//	int i=0;
//	double u_sol, uTs;	
/*
   PXYZ pxy[12];

   double u[12];
   double uT[12];
   double uu[12];
   double uv[12];

   double uu_sol,uv_sol, ux, uy, vx, vy;
   const double c = CShapePara; */

//	int im, iy, ix, iy1, ix1;

//	int ils,jls;  

//	for(i=0; i<NPOINTCIRCLE; i++)
//	{
//		im=0;
//		iy = Company[i].iy;
//		ix = Company[i].ix;
//		iy1 = Company[i+1].iy;
//		ix1 = Company[i+1].ix;
/*
   for(ils=0; ils<3; ils++)
   {
   for(jls=0; jls<3; jls++)
   {					
   if(bodygrid[iy-1+ils][ix-1+jls].flag % 2 == 0)
   {
   pxy[im].x = bodygrid[iy-1+ils][ix-1+jls].xc1;
   pxy[im].y = bodygrid[iy-1+ils][ix-1+jls].yc1;

   u[im] = bodygrid[iy-1+ils][ix-1+jls].jbnp.p;
   uT[im] = bodygrid[iy-1+ils][ix-1+jls].jbnp.T;
   uu[im] = bodygrid[iy-1+ils][ix-1+jls].jbnp.u;
   uv[im] = bodygrid[iy-1+ils][ix-1+jls].jbnp.v;
   im++;
   }
   }
   }
//============
uu[im]=uv[im]=0.0;
pxy[im]=Circle[i];
///==============================
rbf_mq(im, pxy, Circle[i], u, c, u_sol);
rbf_mq(im, pxy, Circle[i], uT, c, uTs);
*/
//		u_sol=0.5*(bodygrid[iy][ix].jbnp.p+bodygrid[iy1][ix1].jbnp.p);  //直接使用Interface Cell中心的值
//		uTs = 0.5*(bodygrid[iy][ix].jbnp.T+bodygrid[iy1][ix1].jbnp.T); //直接使用Interface Cell中心的值

//		Company[i].pre = u_sol;	
//		Company[i].cp=(u_sol-initP)/(0.5*initQ);
//		Company[i].T=uTs;

/*		rbf_mqUxUy(im+1, pxy, Circle[i], uu, c, uu_sol, ux, uy);
        rbf_mqUxUy(im+1, pxy, Circle[i], uv, c, uv_sol, vx, vy); 

        double mu, tx, ty, nx, ny;
        mu = 1.458e-6/FreeMu*pow(FreeT*uTs, 1.5)/(FreeT*uTs+110.4);

        u_sol = sqrt(Company[i].normS.x*Company[i].normS.x
        +Company[i].normS.y*Company[i].normS.y);
        nx = Company[i].normS.x / u_sol;
        ny = Company[i].normS.y / u_sol;
        tx = ny;
        ty = -nx;
        Company[i].cf = 2*mu/Rel * (ux*tx*nx+vx*ty*nx + uy*tx*ny + vy * ty *ny);*/
//	}
//	Company[NPOINTCIRCLE] = Company[0];
//}
/*
   void computeNormalArea()
   {
   int i=0;

   for(i=0; i<NPOINTCIRCLE; i++)
   {
   Company[i].normS.x = Circle[i+1].y-Circle[i].y;
   Company[i].normS.y = Circle[i].x-Circle[i+1].x;
   }
//++++++

//++++++
}

extern double CL,Cd;

void computeCL_Cd()
{
int i;

double sum_fx=0.0;
double sum_fy=0.0;

for(i=0;i<NPOINTCIRCLE;i++)
{
sum_fx += Company[i].normS.x  * Company[i].pre;
sum_fy += Company[i].normS.y  * Company[i].pre;
}

CL = 2.0*sum_fy;
Cd = 2.0*sum_fx;
}


*/





//np:内部插值节点数目
//nb:Neumann边界条件节点数目
//pxy:np个内部节点坐标
//pxynb:nb个Neumann边界条件节点坐标
//xy:待求的点坐标
//u[np]:内部节点的值
//unb[nb]:Neumann边界条件值
//nml[nb]:nb个法向方向
//u_sol:插值结果
/*
   void rbf_mq_nb( int np, int nb, PXYZ *pxy, PXYZ *pxynb, PXYZ &xy, PXYZ *nml,
   double *u, double *unb, double c, double &u_sol);

   void inv_dis_nb( int np, int nb, PXYZ *pxy, PXYZ *pxynb, PXYZ &xy, PXYZ *nml,
   double *u, double *unb, double c, double &u_sol);
//np: 插值节点点数
//pxy [np] 0-------np-1
//u[np]:np个插值节点上的函数值
//c:coefficient
//u_sol插值结果
void rbf_mq( int np, PXYZ *pxy, PXYZ &xy, double *u, double c, double &u_sol);

void shblToJbbl(SHBL *ush, JBBL *jbu);
void jbblToShbl(JBBL *jbu, SHBL *ush);

void inv_dis_nb( int np, int nb, PXYZ *pxy, PXYZ *pxynb, PXYZ &xy, PXYZ *nml,
double *u, double *unb, double c, double &u_sol)
//come from JCP 225, 2007, 2098-2117
{
int i=0;
double sum_af=0.0, sum_bt=0.0, sum0;
double *paf=NULL; 
double *pbt=NULL;

paf = new double [np];
pbt = new double [nb]; 

for(i=0; i<np; i++)
{
paf[i]=1.0 / distance2p(xy, pxy[i]);
sum_af += paf[i];
}

for(i=0; i<nb; i++)
{
pbt[i]=1.0 / distance2p(xy, pxynb[i]);
sum_bt += pbt[i];
}

sum0 = sum_af + sum_bt;

sum_af = 0.0;

for(i=0; i < np; i++)
{
sum_af += paf[i] / sum0 * u[i];
}

u_sol = sum_af/(1.0 - sum_bt/sum0);

delete [] paf;
delete [] pbt;
}


void getInterfaceCell(Node *pnode)
{
const int nmax1=5;
const int nmax2=45;

const double c=CShapePara;

PXYZ wpoint[nmax1];
PXYZ nml[nmax1];
//OctCell *grid[nmax2]={NULL};
int num, numnb, np;

int iy, ix;

//	findInterpolatePoint(iy, ix, grid, num, wpoint, numnb, nml);
num=pnode->num;
numnb=pnode->numnb;

//dirichelet boundary condition  
np = num + numnb;

PXYZ *pxy=new PXYZ[num+numnb];
double *u0=new double [num+numnb];
double *u1=new double [num+numnb];

double usol_p, usol_u, usol_v, usol_T, rho, usol_k, usol_o;
double dd=0.0, d1=0.0, omgwall=0.0, mu=0.0;

PXYZ xy;
int i;

for(i=0; i<np; i++)
{
    if(i<num)
    {
        iy=pnode->pigrid[i][0];
        ix=pnode->pigrid[i][1];

        pxy[i].x=bodygrid[iy][ix].xc1;
        pxy[i].y=bodygrid[iy][ix].yc1;
        u0[i]=bodygrid[iy][ix].jbnp.u;
        u1[i]=bodygrid[iy][ix].jbnp.v;
        mu += bodygrid[iy][ix].mu;
    }
    else
    {
        pxy[i].x=pnode->pwnb[i-num].x;
        pxy[i].y=pnode->pwnb[i-num].y;

        wpoint[i-num]=pnode->pwnb[i-num];

        u0[i]=0.0;
        u1[i]=0.0;
    }
}

mu = mu/num;

xy.x = bodygrid[pnode->iy][pnode->ix].xc1;
xy.y = bodygrid[pnode->iy][pnode->ix].yc1;

rbf_mq( np, pxy, xy, u0, c, usol_u);	
rbf_mq( np, pxy, xy, u1, c, usol_v);

for(i=0; i<np; i++)
{
    if(i<num)
    {
        iy=pnode->pigrid[i][0];
        ix=pnode->pigrid[i][1];

        u0[i]=bodygrid[iy][ix].jbnp.k;
    }
    else
    {
        u0[i]=0.0;
    }
}

rbf_mq( np, pxy, xy, u0, c, usol_k);

//void rbf_mq_nb( int np, int nb, PXYZ *pxy, PXYZ *pxynb, PXYZ &xy, PXYZ *nml,
//	   double *u, double *unb, double c, double &u_sol);

double *unb0=new double [numnb];

for(i=0; i<num; i++)
{
    iy=pnode->pigrid[i][0];
    ix=pnode->pigrid[i][1];

    u0[i]=bodygrid[iy][ix].jbnp.p;
    u1[i]=bodygrid[iy][ix].jbnp.T;
}

for(i=0; i < numnb; i++)
{
    unb0[i]=0.0;
}

#ifdef RBF_NEUMANN
rbf_mq_nb( num, numnb, pxy, wpoint, xy, nml,
    u0, unb0, c, usol_p);

rbf_mq_nb( num, numnb, pxy, wpoint, xy, nml,
    u1, unb0, c, usol_T);
#elif defined INVERSE_DISTANCE_NEUMANN
inv_dis_nb( num, numnb, pxy, wpoint, xy, nml,
    u0, unb0, c, usol_p);

inv_dis_nb( num, numnb, pxy, wpoint, xy, nml,
    u1, unb0, c, usol_T);
#endif 

rho = usol_p/usol_T/Rstar;

dd=0.0;
for(i=0; i< numnb; i++)
{
    dd += sqrt((wpoint[i].x-xy.x) * (wpoint[i].x-xy.x) 
        + (wpoint[i].y-xy.y) * (wpoint[i].y-xy.y));
}

d1 = dd/(numnb+1e-18);

omgwall=60.0*FreeMu*mu/rho/FreeD/0.075/FreeOmg/d1/d1/Length/Length;
for(i=0; i<np; i++)
{
    if(i<num)
    {
        iy=pnode->pigrid[i][0];
        ix=pnode->pigrid[i][1];
        u0[i]=bodygrid[iy][ix].jbnp.omg;
    }
    else
    {
        u0[i]=omgwall;
    }
}
rbf_mq( np, pxy, xy, u0, c, usol_o);

delete [] u0;
delete [] u1;
delete [] unb0;
delete [] pxy;

JBBL jbu;
SHBL shu;

jbu.p=usol_p;
jbu.u=usol_u;
jbu.v=usol_v;
jbu.T=usol_T;
jbu.q=rho;
jbu.k=usol_k;
jbu.omg=usol_o;

bodygrid[pnode->iy][pnode->ix].jbnp=jbu;
bodygrid[pnode->iy][pnode->ix].jbtmpr=jbu;

jbblToShbl(&jbu, &shu);

bodygrid[pnode->iy][pnode->ix].shbl1np = shu;
bodygrid[pnode->iy][pnode->ix].sh1tmpr = shu;
bodygrid[pnode->iy][pnode->ix].sh1tmpr0 = shu;
}

void getInternalSolidCell(int iy,int ix);


void wall_treatment()
{

    Node *current;
    current = HeadFaceCell;
    while(current != NULL)
    {
        getInterfaceCell(current);

        current = current->next;
    }
}

*/



