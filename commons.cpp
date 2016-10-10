
#include"non_uniform_grid.h"

OctCell *bodygrid;

//OctCell *refcell[Nr+1][NneedRefineGrid];   //--记录需要加密的网格单元--
vector<OctCell *> refcell[Nr+1];
//OctCell *solrefcell[NneedRefineGrid];
vector<OctCell *>solrefcell;//solution AMR cell

vector<OctCell *>solcoacell;//solution coase cell

PXYZ MultiCircle[N_BODY+1][500]; //points at solid wall
int NWallPts[N_BODY+3]={NPOINTCIRCLE, NPOINTCIRCLE1,
    NPOINTCIRCLE2};

string WallPtsFile[N_BODY+3]={ WALL_POINT_FILE, WALL_POINT_FILE1,
 WALL_POINT_FILE2};
//PXYZ Circle[NPOINTCIRCLE+1]; //points at solid wall
//PXYZ Circle1[NPOINTCIRCLE1+1]; //points at solid wall
//PXYZ Circle2[NPOINTCIRCLE2+1]; //points at solid wall
PXYZ ext_wall[NPT_EXT_WALL+1]; //points at external computational domain

Node *HeadListAllGrid=NULL; //the head of all cell 

BndryNode *HeadFaceCell=NULL;//the head of boudnary internal cell 
vector<OctCell*> ExtWallBndry;
vector<OctCell*> InletBndry;
vector<OctCell*> OutletBndry;

vector<Face> faces; //临时面

vector<Face> faces_comp;  //需要计算的面

map<OctCell *, cellEdgeBool> cellEbools; //用来找边的临时map

map<OctCell *, cellPointIndex> cellPoints;  //cell由那几个点组成序号为points4out的下标+1

vector<pointInCell> points4out; //输出时候用,给出点的坐标

double dt; //global time step
double timesum=0.0;

set<OctCell*>mvalue_set; //可能的多值cell集合
PXYZ MBsharp_point[N_BODY+1][8];
//多值点个数不止一个时候不要忘记solveAMR.cpp中还有个确定多值点的函数不要少调用了
Vec2D ExtPnt0(-3.0,2.0);
Vec2D ExtPnt1(-3.0,-2.0);
Vec2D ExtPnt2(5.0,-2.0);
Vec2D ExtPnt3(5.0,2.0);
double CL,Cd;

double value_kxrcf_amr=1.0;
string str0("OUT2//");
JBBL State_l;
JBBL State_r;

set<OctCell*>ColCoa;
double sumVorM;
double epsPos=1.0e-13;
double max_cv_x;
double max_cv_y;
double CFL=0.10;
