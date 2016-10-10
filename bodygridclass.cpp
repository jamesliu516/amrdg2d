
#include "non_uniform_grid.h"
// defaut construction function

node::node() 
{
    cell=NULL;
    next=NULL;
    prev=NULL;
    flg=0;
#ifdef RES_SMOOTH
    n_nb=0;
#endif
}
bodygridclass::bodygridclass()
{
    //   pmvx=NULL;
    //   pmvy=NULL;

    flag=0;
    exflag=-1;
    level=0;//网格所在的层数 
    reflag=0;//是否加密;0-未加密,1-已加密
    coflag=false;//判断网格是否粗化；0-未粗化,1-已粗化
    NOparent=0u;//初始网格的序号
    NOchildren=0;//该网格单元在子网格中的编号
    parent=NULL; 
    bAvg=false;
    tcur=0.0; 
    tdiv=0.0;
    pio=0.0;
    trb=false;
    trbamr=false;
    kxq=kxrcf=0.0;
    maxqp=inlen=0.0;
    mvflg=false;
    vorM=0.0;
    //	struct bodygrid *children[4];

//    p_node=NULL;
    for(int m=0; m<4; ++m)  children[m]=NULL;

    for(int i=0; i<4; i++) { 
        for(int j=0; j<6; j++){
#ifdef RES_SMOOTH
            resd0[i][j]=resd1[i][j]=0.0;
#endif
            residual[i][j]=dof[i][j]=dof0[i][j]=0.0;
        }
    }
}
//带参数构造函数
bodygridclass::bodygridclass(double hhx, double hhy)
{
    //   pmvx=NULL;
    //   pmvy=NULL;

    flag=0;
    exflag=-1;
    level=0;//网格所在的层数 
    reflag=0;//是否加密;0-未加密,1-已加密
    coflag=false;//判断网格是否粗化；0-未粗化,1-已粗化
    NOparent=0u;//初始网格的序号
    NOchildren=0;//该网格单元在子网格中的编号
    parent=NULL; 
    bAvg=false;
    tcur=0.0; 
    tdiv=0.0;
    pio=0.0;
    trb=false;
    trbamr=false;
    kxq=kxrcf=0.0;
    maxqp=inlen=0.0;
    mvflg=false;
    vorM=0.0;
    //	struct bodygrid *children[4];

   // p_node=NULL;
    for(int m=0; m<4; ++m)  children[m]=NULL;

    for(int i=0; i<4; i++) { 
        for(int j=0; j<6; j++){
#ifdef RES_SMOOTH
            resd0[i][j]=resd1[i][j]=0.0;
#endif
            residual[i][j]=dof[i][j]=dof0[i][j]=0.0;
        }
    }
    set_invM(hhx,hhy);
}

Face::Face()
{
    parent=NULL; //left cell
    neighbor=NULL; //right cell
    //parent is the id of the cell owning this face
    // A cell owns a face if the face normal is pointing outwards from it
    // The other cell is called the neighbor
    //  Vec2D nml;
    nml[0]=0.0;
    nml[1]=0.0; 

    bHang=false;
    for(int ii=0; ii<3; ++ii){
        gs_pt[ii][0]=gs_pt[ii][1]=0.0;
        for(int jj=0; jj<4; ++jj)
            faceflux[ii][jj]=0.0; //face flux at 3 gauss points
    }         
    area=0.0; //the area of the face in 3d or lenth in 2d
}

cellEdgeBool::cellEdgeBool()
{
    for(int i=0; i<4; ++i)
        bEdge[i]=false;
}

cellPointIndex::cellPointIndex()
{
    for(int i=0; i<4; ++i)
        iPoint[i]=0;   
}


pointInCell::pointInCell()
{
    nmc=0;
    for(int i=0;i<4;++i) mcell[i]=NULL;
}


