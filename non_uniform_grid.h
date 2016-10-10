#ifndef NS_AMR_H
#define NS_AMR_H  

#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>
#include<vector>
#include<set>
#include<list>
#include<map>
#include<utility>
#include<algorithm>

#include"vec2d.h"
#include"vec3d.h"
#include"vec4d.h"

using namespace std;

#define DEBUG

#define WALL_BNDRY
//#define PERIODIC_BNDRY          //periodic boundary
//#define DOUBLE_MA_BNDRY // double mach reflection problem
//#define RIEMANN2_BNDRY

//#define SHBL_LIM
#define TZBL_LIM
const int HLLCFLUX=0; //1 hllc,0 global Lax-Friedrichs flux

//#define DUICHEN  
//#define ILWILW

//#define RES_SMOOTH  //residual smoothing
//#define LOCAL_TIME //局部时间步加速
//#define FIXED_DT  //fixed time step
const double fixed_dt=1e-5;
#define FILE_PROX  "dg"

//const double CFL= 0.10444;  //0.4

const double t_print=0.245;//指定时间处的输出

const int Limiter01=1; 
const double M4TVBM =1.010;
#define KXRCF_LIM
const int PosLimiter=1;// 1 with positive preserving, 0 no 
const int NLBpt=3; //3 Gauss LB pts, 4 gass LB pts, others is illegal
//const double lim_muscl=1.0;  //=0.5 or 1.0

const int nDOF=6; //6: for p2, 3 for p1 其他数据出错
//#define NORMAL_DIR_MOD  //边界法向修正暂时不能用于有奇点的外形
const int qvlvModify=1; //0直接使用对称边界，1使用曲率修正

const bool b_ARS=true;//do use riemann soleve to set wall bound
const int qvlvRie=0;
const int mv_qvlvModify=0;//多值点附近只是对称边界条件
const double xLeft= -3.3;
const double xRight= 5.3;
const double yLow= -2.2;
const double yUpper= 2.2; 

const double xInlet=-3.0;
const double xOutlet=5.0;
const double PI = 3.14159265358979;

const int Nx=123;  //total cell number in x direction 其中左边两个以及右边两个用于边界设置
const int Ny=65;//详细说明也可以参考formlist.cpp
const int max_time_steps=80000;
const bool isOutVorM=false;
//**********************************************************
const bool SwAMR=true; //是否做解自适应 
const int Nr= 3;  // 初始网格初始几何加密Nr次
const int Nar=4; //解自适应次数
const int MaxAMR=4;
#define NAR_MAXAMR 0   //1 use Nar to assure the number of AMR
//0 use MaxAMR to make sure the number of AMR
//#define KXRCF_AMR
const int AMRP=4; //solution AMR refine control:
const short COAR=4;//amr indicator choose for coase 
//0:cur div, 1:div, 2:curl, 3:pio, others, 4:other combination
//AMRP=5 KXRCF with div
//对于3:other 自适应标度
#define ADAPTIVE_INDICATOR 0 //0:den,1:press, 2:entropy
const double coarse_coe=0.3; //amr coarse coefficient can be changed
const double coarse_coe1=0.4; //amr coarse coefficient can be changed
const double refine_coe=1.2;
const double refine_coe1=1.60; //次要加密设置参数
const int localRefine=0; //一般不用这个除非就是直接使用
//局部加密Nr=1,2,3,4,5其他数值要重新写外部区域加密的程序??

const int adp_begin=7;//解自适应从每个坐标方向第adp_begin+1个开始,大一点可以防止远场边界设置会出现内点不是叶子节点的问题
//*************************************************************
#define GRID_FILE  "IN//naca0012Grid2.in"
//#define WALL_POINT_FILE  ("IN//RAE2822.dat")
//#define WALL_POINT_FILE  ("IN//halfCircle")
#define WALL_POINT_FILE  ("IN//tri_Pro")
//#define WALL_POINT_FILE  "IN//wallpoint1.dat"
#define WALL_POINT_FILE1  "IN//wallpoint2.dat"
#define WALL_POINT_FILE2  "IN//naca0012wall.in"
//#define WALL_POINT_FILE  "IN//scramjet_point"
#define ExtWALL_POINT_FILE  "IN//quad01_point"
#define IS_TUBE 1 //1 yes 0 no
#define N_BODY 1 //!=0,  1:  one body 2:  two bodies, 3: three bodies and so on. 
#define HAVE_MV 1// 1 with multi value points, 0 no
//多值点标号与个数
const int n_sharp_point[N_BODY+3][8]={{0,100,200},{78}};
const int MVptNum[N_BODY+3]={3};
//#define WALL_POINT_FILE  ("IN//naca0012bd.dat")
/*+++++++++++++++++++++++++++++++++++++++++++++++++
  circleGrid.in : Nx=298+2, Ny=338+2, NPOINT=102641, NCELL=102000
  circleGrid1.in : Nx=445+2, Ny=255+2, NPOINT=115584, NCELL=114879
  circleGrid2.in : Nx=535+2, Ny=325+2, NPOINT=176464, NCELL=175599
  naca0012Grid.in : Nx=285+2, Ny=296+2 , NPOINT= 86112, NCELL=85526
  naca0012Grid2.in : Nx=285+2, Ny=296+2 , NPOINT= 86112, NCELL=85526
  naca0012Grid3.in : Nx = 535+2 Ny = 716+2, NPOINT= 386822, NCELL=385566
  +++++++++++++++++++++++++++++++++++++++++++++++++*/
const int NPOINT= 86112;
const int NCELL= 85526;

const int NPOINTCIRCLE = 300;  // the point at the wall for circel =360,
const int NPOINTCIRCLE1 = 154;  // the point at the wall for circel =360,
const int NPOINTCIRCLE2 = 154;  // the point at the wall for circel =360,
// for naca0012 = 154, no point is same. for RAE2822:118
const int NPT_EXT_WALL=312;//points at external compuatianial domain 

const double FreeMa =10.0; //来流马赫数

const double AOA=0.0;   //incidence or angle of attack单位度
//此时速度的方向改变，就是体轴系
//物体旋转角度 单位度数
const double rotate_deg=AOA;//换成风轴系

//dimensionless constant 来流用国际单位
const double FreeD=1.225;  // free density use d,q to denote
const double FreeT=288.15;     //free temperature

const double GAMMA=1.4;
const double GAM11=GAMMA-1.0;

const double Rconst = 287.04;
const double Rstar=1.0/GAMMA; //使用音速作为速度的无量纲方式，可以参考阎超的书
const double FreeP=FreeD*Rconst*FreeT; //来流用国际单位
const double FreeV=FreeMa*340.29;        //free velocity来流用国际单位
const double FreeA=sqrt(GAMMA*FreeP/FreeD);

const double initP = FreeP/(FreeD*FreeA*FreeA);   //---无量纲标准大气压 音速a=340.29-
const double initQ = 1.0;   //---无量纲标准大气密度---
const double initT=initP/initQ/Rstar;

const double uFree=FreeMa;  //无量纲速度
const double vFree=0.0;
//const double uFree=FreeMa*cos(PI/180.0 * AOA);  //无量纲速度
//const double vFree=FreeMa*sin(PI/180.0 * AOA);

const double hx=(xRight-xLeft) / (Nx-4.0);  //等距网格//这里用来表示等间距网格的情形
const double hy=(yUpper-yLow) / (Ny-4.0);

const double RADII= 0.300;  //无量纲后的量  圆的半径
const double x_o=0.8;
const double y_o=0.0;//圆心坐标

const double Ellipse_x=-0.02; //局部加密的位置指定开始加密几层,局部指定位置加密时需要使得Nr=0
const double Ellipse_x1=1.02; //局部加密的位置指定开始加密几层,局部指定位置加密时需要使得Nr=0
const double Ellipse_y=-0.58;
const double Ellipse_y1=0.58;

const double ERRS = 1.0e-8;

const double gspt[3]={-0.774596669241483, 0.00000000000, 0.774596669241483};
const double wei[3]={0.555555555555556, 0.888888888888889, 0.555555555555556};  

const double gsLB3pt[3]={-1.0,0.0,1.0};
const double weiLB3[3]={1.0/3.0, 4.0/3.0, 1.0/3.0};

const double gsLB4pt[4]={-1.0, -0.447213595499958,
    0.447213595499958, 1.0};
const double weiLB4[4]={1.0/6.0, 5.0/6.0, 5.0/6.0, 1.0/6.0};

//const double rk_coe[4]={};
//#define VTK_OUTPUT
const double RXMIN = 1.5;
const double RXMAX =7.5;

const double RYMIN = 0.15; //-0.065 for airfoil
const double RYMAX =0.7;   //0.065 for airfoil
#define TECPLOT_OUTPUT

#define max2(a,b) ((a)>(b)?(a):(b))
#define max3(a,b,c) max2(max2(a,b),c)
#define min2(a,b) ((a)<(b)?(a):(b))
#define min3(a,b,c) min2(min2(a,b),c)
//#define diffU(a,b) ((a)-(b))

//目前都用不着了
#define LEN sizeof(OctCell)
#define NpointForOutput 800000
#define NcellForOutput 800000
#define NneedRefineGrid 800000
#define NMAX 40   //用于二维方程组求解的程序
const int reserve_num=400000;
typedef struct conservedvariables{
    double q;
    double qu;
    double qv;
    double te;   //---单位体积总能E=q*(e+0.5*V*V), e=p/(gama-1)q---

    conservedvariables() {q=qu=qv=te=0.0;}

}SHBL;

typedef struct privariables{
    //	double T;
    double q;
    double u;
    double v;
    double p; //--total specific enthalpy 总焓H=h+0.5*V*V,h=e+p/q---

    privariables(){u=v=q=p=0.0;}
}JBBL;

typedef class pointxy{
public:
    double x;
    double y;
    pointxy(){x=y=0.0;}
}PXYZ;

class bodygridclass;

typedef class wallp
{
public:
    Vec2D normS;
    double pre;
    //   int iy;
    // int ix;
    bodygridclass *incell;
    double cp;
    double cf;
    double T;
    wallp()
    {
        pre=0.0;
        incell=NULL;
        cp=cf=0.0;
        T=0.0;
    }
}WALLP;


//typedef struct flux{
//	double fq;
//	double fqu;
//	double fqv;
//	double fte;
//}Flux;

class node;

class BndryNode;

typedef class bodygridclass
{ 
public:
    double invM[6]; //local mass matrix
    double dof[4][6];
    double dof0[4][6];
    short flag; //网格的类型：0-fluid cell, -1-solid cell, 2-interface cell 
    short exflag;//用于内流固壁边界条件的设置0-fluid, -1 solid, 2 interface cell 
    //flag:当网格中心在固体内部且为离壁面最近1层的cell设为2，此用于设置边界条件，
    // 中心在固体内部的其他网格设为-1,不参与设置虚拟点，及计算
    //网格中心在流体中设为0  2011-9-8
    //---------------------------------------
//multi value flag 可能的不一定就是
    bool mvflg;
    double Ux[4],Uy[4];//gradient
    short level;//网格所在的层数 
    short level0;//网格在自适应之前的初始层数
    short reflag;//是否加密;0-未加密,1-已加密
    bool coflag;//判断网格是否粗化；0-未粗化,1-已粗化

    unsigned NOparent;//初始网格的序号

    short NOchildren;//该网格单元在子网格中的编号

    bool bAvg;//用于说明是否设置了单元平均
    double vorM; //vortex

    //自适应判据
    bool trb;//trouble cell true or false
    bool trbamr;//trouble cell true or false for AMR
    double kxrcf;
    double kxq;
    double maxqp;
    double inlen;
    double tcur;
    double tdiv;
    double pio;//all other adaptive indicator
//limiter for positive preserving

    bodygridclass *parent; 
    bodygridclass *children[4];
	/*---------------------------------
			         3     2
			         0     1
	              7     6
		     4     5
    ---------------------------------*/

    double xc1,yc1;//网格单元中心的坐标
    double dx,dy; //网格步长

    double dt1;   //local time step

    double residual[4][6];
#ifdef RES_SMOOTH
    double resd0[4][6];
    double resd1[4][6];
#endif

 //   node *p_node;
    //---------------------------------------------------

// BndryNode *pnode;  //首先将那些不是多值点的cell先设置好,边界网格处理时候如果不是多值cell直接取值,这时pnode=NULL
    //要是多值cell在时间推进过程中计算ghost控制体值pnode!=NULL

public:
    bodygridclass();
    bodygridclass(double,double);

    void set_invM() {
        invM[0]=1.0/dx/dy;   
        invM[1]=3.0/dx/dy;
        invM[2]=3.0/dx/dy;
        invM[3]=9.0/dx/dy;
        invM[4]=45.0/4.0/dx/dy;
        invM[5]=45.0/4.0/dx/dy;
    }

    void set_invM(double dsx, double dsy) {
        invM[0]=1.0/dsx/dsy;   
        invM[1]=3.0/dsx/dsy;
        invM[2]=3.0/dsx/dsy;
        invM[3]=9.0/dsx/dsy;
        invM[4]=45.0/4.0/dsx/dsy;
        invM[5]=45.0/4.0/dsx/dsy;
    }

}OctCell;  

typedef OctCell *ptrOctCell;

typedef class node
{
public:
    node();
    short flg;  //flg=6 farfied cell用来设置边界条件的, 
    //flg=2 除了可以自适应的cell以外的一部分流体需要计算的cell
    //flg=0 cell can be adaptive

    OctCell *cell;        
    node *next;
    node *prev;              
#ifdef RES_SMOOTH
    short n_nb;
#endif
/*    ~node()
    { 
        if(cell!=NULL) cell->p_node=NULL;
    }
    */
}Node;  

class BndryNode {
public:
    BndryNode()
    {
        cell=NULL;
        next=NULL;
        prev=NULL;
     //   incell=NULL;
        //    flg=0;    
    }

    //   short flg; //flg=2一般固体内部cell设置一般边界条件, flg=4多值点cell                  
    OctCell *cell; 

    BndryNode *next;
    BndryNode *prev;  

   // PXYZ smtrc_pt;//symmetrical point

    //OctCell *incell; //symmetrical point 所在的cell
    //vector<OctCell*> nbcell; //symmetrical point 邻居插值cell

};

class Face {
public:
    OctCell *parent; 
    OctCell *neighbor; 
    //parent is the id of the cell owning this face
    // A cell owns a face if the face normal is pointing outwards from it
    // The other cell is called the neighbor
    // Vec2D nml;
    double nml[2];
    double gs_pt[3][2];//存放gauss积分点坐标
    double faceflux[3][4]; //face flux at 3 gauss points
    double area; //the area of the face in 3d or lenth in 2d
    bool bHang;
    Face();
};

class cellEdgeBool {
public:
    bool bEdge[4];
    cellEdgeBool();

};    

//cell由哪四个点组成
class cellPointIndex {
public:
    int iPoint[4];  //initial初始化为 0 表示点的序号从1开始的
    cellPointIndex();
};

//点for output
class pointInCell {
public:
    Vec2D pt;
    OctCell *mcell[4];//包含这个点的cell,最多4个,输出时用的
    short nmc;
    pointInCell();
};

inline void shblToJbbl(SHBL *ush, JBBL *jbu)
{
    jbu->q=ush->q;
    jbu->u=ush->qu/ush->q;
    jbu->v=ush->qv/ush->q;

    jbu->p = (ush->te - 0.5 * ush->q * (jbu->u * jbu->u
            + jbu->v * jbu->v)) * GAM11;

}

inline void shbl2jbbl(double ush[],double ujb[])
{
    ujb[0]=ush[0];
    ujb[1]=ush[1]/ush[0];
    ujb[2]=ush[2]/ush[0];
    ujb[3]=(ush[3]-0.5*ush[0]*(ujb[1]*ujb[1]+ujb[2]*ujb[2]))*GAM11;
}

inline void jbblToShbl(JBBL *jbu, SHBL *ush)
{
    ush->q=jbu->q;
    ush->qu=jbu->u * jbu->q;
    ush->qv=jbu->v * jbu->q;

    ush->te=jbu->p/GAM11+ 0.5 * ush->q * ( jbu->u * jbu->u
        + jbu->v * jbu->v );

}
inline void jbbl2shbl(double ujb[],double ush[])
{
    ush[0]=ujb[0];
    ush[1]=ush[0]*ujb[1];
    ush[2]=ujb[0]*ujb[2];
    ush[3]=ujb[3]/GAM11+0.5*ujb[0]*(ujb[1]*ujb[1]+ujb[2]*ujb[2]);
}

// base function
inline double  phix(double x, double xc1,  double dx)
{
    return (2.0*(x-xc1)/dx);
}

inline double psiy(double y, double yc1, double dy)
{
    return (2.0*(y-yc1)/dy);
} 

inline double phix_psiy(double x,double y, double xc1, double yc1, double dx, double dy)
{
    return (4.0*(y-yc1)*(x-xc1)/(dx*dy));
}   

inline double phix2m(double x, double xc1, double dx)
{
    return (4.0*(x-xc1)*(x-xc1)/(dx*dx)-1.0/3.0);
}

inline double psiy2m(double y,  double yc1, double dy)
{
    return (4.0*(y-yc1)*(y-yc1)/(dy*dy)-1.0/3.0);
}

inline double distance2p(const PXYZ &p1, const PXYZ &p2)
{
    return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

#endif

