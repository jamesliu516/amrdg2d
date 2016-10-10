
#include"non_uniform_grid.h"
#include"findneighbor.h"

class Facetmp{
public:
    OctCell *parent;
    OctCell *neighbor;

    Facetmp(){parent=neighbor=NULL;}

    bool operator==(const Facetmp&);
};

bool Facetmp:: operator==(const Facetmp& right1)
{
    return ((parent==right1.parent && neighbor==right1.neighbor) || 
        (parent==right1.neighbor && neighbor==right1.parent) );
}

extern map<OctCell *, cellEdgeBool> cellEbools;

extern vector<Face> faces; //all face
extern vector<Face> faces_comp;  //faces will be computed 
extern Node *HeadListAllGrid;

//     2
//    ______
//   |     | 
// 3 |     | 1
//   ------
//     0

//OctCell *NorthNeighbor(OctCell *pp);
//OctCell *WestNeighbor(OctCell *pp);
//OctCell *EastNeighbor(OctCell *pp);
//OctCell *SouthNeighbor(OctCell *pp);


void face_vec()
{  
    //map<OctCell *, cellEdgeBool>::size_type cnt0;

    OctCell *pcell0;
    OctCell *pcell1;   
    Face fff0;

    Node *current;
    current = HeadListAllGrid;   

    //   int cnt=0;

    faces_comp.clear();
    while(current != NULL) {
        //-------
        if(current->flg<=2) {
            pcell0=current->cell;             
            //south neighbor
            pcell1=SouthNeighbor(pcell0);           
            if(pcell1->reflag==0) {
                if(pcell0->level==pcell1->level) {
                    if(cellEbools[pcell0].bEdge[0]==false
                        && cellEbools[pcell1].bEdge[2]==false){ 
                        cellEbools[pcell0].bEdge[0]=true;
                        cellEbools[pcell1].bEdge[2]=true;
                        fff0.parent=pcell0;
                        fff0.neighbor=pcell1;
                        fff0.nml[0]=0.0;
                        fff0.nml[1]=-1.0;

                        fff0.area=pcell0->dx;

                        fff0.gs_pt[0][1]=fff0.gs_pt[1][1]=fff0.gs_pt[2][1]=pcell0->yc1-0.5*pcell0->dy;
                        for(int ii=0;ii<3;++ii) fff0.gs_pt[ii][0]=gspt[ii]*0.5*pcell0->dx+pcell0->xc1;

                        faces.push_back(fff0);
                    }
                }
                else if(pcell0->level!=pcell1->level){
                    if(cellEbools[pcell0].bEdge[0]==false) {
                        cellEbools[pcell0].bEdge[0]=true;
                        fff0.parent=pcell0;
                        fff0.neighbor=pcell1;
                        fff0.nml[0]=0.0;
                        fff0.nml[1]=-1.0;

                        fff0.bHang=true;
                        fff0.area=pcell0->dx;

                        fff0.gs_pt[0][1]=fff0.gs_pt[1][1]=fff0.gs_pt[2][1]=pcell0->yc1-0.5*pcell0->dy;
                        for(int ii=0;ii<3;++ii) fff0.gs_pt[ii][0]=gspt[ii]*0.5*pcell0->dx+pcell0->xc1;

                        faces.push_back(fff0);                       
                    } 
                }
            }  


            //east neighbor
            pcell1=EastNeighbor(pcell0);           
            if(pcell1->reflag==0) {
                if(pcell0->level==pcell1->level) {
                    if(cellEbools[pcell0].bEdge[1]==false && cellEbools[pcell1].bEdge[3]==false){ 
                        cellEbools[pcell0].bEdge[1]=true;
                        cellEbools[pcell1].bEdge[3]=true;
                        fff0.parent=pcell0;
                        fff0.neighbor=pcell1;
                        fff0.nml[0]=1.0;
                        fff0.nml[1]=0.0;

                        fff0.area=pcell0->dy;
                        fff0.gs_pt[0][0]=fff0.gs_pt[1][0]=fff0.gs_pt[2][0]=pcell0->xc1+0.5*pcell0->dx;
                        for(int ii=0;ii<3;++ii) fff0.gs_pt[ii][1]=gspt[ii]*0.5*pcell0->dy+pcell0->yc1;  

                        faces.push_back(fff0);
                    }
                }
                else if(pcell0->level!=pcell1->level){
                    if(cellEbools[pcell0].bEdge[1]==false) {
                        cellEbools[pcell0].bEdge[1]=true;
                        fff0.parent=pcell0;
                        fff0.neighbor=pcell1;
                        fff0.nml[0]=1.0;
                        fff0.nml[1]=0.0;
                        fff0.area=pcell0->dy;
                        fff0.bHang=true;

                        fff0.gs_pt[0][0]=fff0.gs_pt[1][0]=fff0.gs_pt[2][0]=pcell0->xc1+0.5*pcell0->dx;
                        for(int ii=0;ii<3;++ii) fff0.gs_pt[ii][1]=gspt[ii]*0.5*pcell0->dy+pcell0->yc1;  

                        faces.push_back(fff0);                       
                    } 
                }
            }            

            //north neighbor
            pcell1=NorthNeighbor(pcell0);           
            if(pcell1->reflag==0) {
                if(pcell0->level==pcell1->level) {
                    if(cellEbools[pcell0].bEdge[2]==false && cellEbools[pcell1].bEdge[0]==false){ 
                        cellEbools[pcell0].bEdge[2]=true;
                        cellEbools[pcell1].bEdge[0]=true;
                        fff0.parent=pcell0;
                        fff0.neighbor=pcell1;
                        fff0.nml[0]=0.0;
                        fff0.nml[1]=1.0;
                        fff0.area=pcell0->dx;                     
                        fff0.gs_pt[0][1]=fff0.gs_pt[1][1]=fff0.gs_pt[2][1]=pcell0->yc1+0.5*pcell0->dy;
                        for(int ii=0;ii<3;++ii) fff0.gs_pt[ii][0]=gspt[ii]*0.5*pcell0->dx+pcell0->xc1;                     

                        faces.push_back(fff0);
                    }
                }
                else if(pcell0->level!=pcell1->level){
                    if(cellEbools[pcell0].bEdge[2]==false) {
                        cellEbools[pcell0].bEdge[2]=true;
                        fff0.parent=pcell0;
                        fff0.neighbor=pcell1;
                        fff0.nml[0]=0.0;
                        fff0.nml[1]=1.0;
                        fff0.bHang=true;
                        fff0.area=pcell0->dx;
                        fff0.gs_pt[0][1]=fff0.gs_pt[1][1]=fff0.gs_pt[2][1]=pcell0->yc1+0.5*pcell0->dy;
                        for(int ii=0;ii<3;++ii) fff0.gs_pt[ii][0]=gspt[ii]*0.5*pcell0->dx+pcell0->xc1;   

                        faces.push_back(fff0);                       
                    } 
                }
            }                  

            //west neighbor
            pcell1=WestNeighbor(pcell0);           
            if(pcell1->reflag==0) {
                if(pcell0->level==pcell1->level) {
                    if(cellEbools[pcell0].bEdge[3]==false && cellEbools[pcell1].bEdge[1]==false){ 
                        cellEbools[pcell0].bEdge[3]=true;
                        cellEbools[pcell1].bEdge[1]=true;
                        fff0.parent=pcell0;
                        fff0.neighbor=pcell1;
                        fff0.nml[0]=-1.0;
                        fff0.nml[1]=0.0;
                        fff0.area=pcell0->dy;
                        fff0.gs_pt[0][0]=fff0.gs_pt[1][0]=fff0.gs_pt[2][0]=pcell0->xc1-0.5*pcell0->dx;
                        for(int ii=0;ii<3;++ii) fff0.gs_pt[ii][1]=gspt[ii]*0.5*pcell0->dy+pcell0->yc1; 

                        faces.push_back(fff0);
                    }
                }
                else if(pcell0->level!=pcell1->level){
                    if(cellEbools[pcell0].bEdge[3]==false) {
                        cellEbools[pcell0].bEdge[3]=true;
                        fff0.parent=pcell0;
                        fff0.neighbor=pcell1;
                        fff0.nml[0]=-1.0;
                        fff0.nml[1]=0.0;
                        fff0.bHang=true;
                        fff0.area=pcell0->dy;                       
                        fff0.gs_pt[0][0]=fff0.gs_pt[1][0]=fff0.gs_pt[2][0]=pcell0->xc1-0.5*pcell0->dx;
                        for(int ii=0;ii<3;++ii) fff0.gs_pt[ii][1]=gspt[ii]*0.5*pcell0->dy+pcell0->yc1; 

                        faces.push_back(fff0);                       
                    } 
                }
            }                
        }
        //========
        //   cout<<cnt<<endl;
        current = current->next;
    }   
    cout<<"Total face: "<<faces.size()<<endl;
    cellEbools.clear();

    for(vector<Face>::size_type isz=0; isz!=faces.size(); ++isz){
        if(!(faces[isz].parent->flag==-1&&faces[isz].neighbor->flag==-1)) {
            faces_comp.push_back(faces[isz]);
        }
    }
    /*
#ifdef RES_SMOOTH
    OctCell *pc0, *pc1;
    for(vector<Face>::size_type isz=0; isz!=faces_comp.size(); ++isz){
            pc0=faces_comp[isz].parent;
            pc1=faces_comp[isz].neighbor;
            if(pc0->flag==0) pc1->p_node->n_nb+=1;
            if(pc1->flag==0) pc0->p_node->n_nb+=1;
    }
#endif
*/
    faces.clear();
    cout<<"Total face will be solved: "<<faces_comp.size()<<endl;      
}



void face_vec1()//±È½ÏÂý
{
    vector<Facetmp> facestmp;

    Facetmp ftmp0;

    OctCell *pcell0;
    OctCell *pcell1;   
    Face fff0;

    Node *current;
    current = HeadListAllGrid;   
    faces.clear();
    //   int cnt=0;
    while(current != NULL) {
        //-------
        if(current->flg<=2) {
            pcell0=ftmp0.parent=current->cell;

            //south neighbor
            pcell1=SouthNeighbor(pcell0);           
            if(pcell1->reflag==0) {
                ftmp0.neighbor=pcell1;                 
                if(find(facestmp.begin(), facestmp.end(), ftmp0 ) == facestmp.end() ) {
                    facestmp.push_back(ftmp0);
                    fff0.parent=pcell0;
                    fff0.neighbor=pcell1;
                    fff0.nml[0]=0.0;
                    fff0.nml[1]=-1.0;
                    faces.push_back(fff0);
                    //         cnt++;
                }
            }

            //east neighbor
            pcell1=EastNeighbor(pcell0);           
            if(pcell1->reflag==0) {
                ftmp0.neighbor=pcell1;                 
                if(find(facestmp.begin(), facestmp.end(), ftmp0 ) == facestmp.end() ) {
                    facestmp.push_back(ftmp0);
                    fff0.parent=pcell0;
                    fff0.neighbor=pcell1;
                    fff0.nml[0]=1.0;
                    fff0.nml[1]=0.0;
                    faces.push_back(fff0);
                    //       cnt++;                    
                }
            }            

            //north neighbor
            pcell1=NorthNeighbor(pcell0);           
            if(pcell1->reflag==0) {
                ftmp0.neighbor=pcell1;                 
                if(find(facestmp.begin(), facestmp.end(), ftmp0 ) == facestmp.end() ) {
                    facestmp.push_back(ftmp0);
                    fff0.parent=pcell0;
                    fff0.neighbor=pcell1;
                    fff0.nml[0]=0.0;
                    fff0.nml[1]=1.0;
                    faces.push_back(fff0);
                    //        cnt++;                     
                }
            }                  

            //west neighbor
            pcell1=WestNeighbor(pcell0);           
            if(pcell1->reflag==0) {
                ftmp0.neighbor=pcell1;                 
                if(find(facestmp.begin(), facestmp.end(), ftmp0 ) == facestmp.end() ) {
                    facestmp.push_back(ftmp0);
                    fff0.parent=pcell0;
                    fff0.neighbor=pcell1;
                    fff0.nml[0]=-1.0;
                    fff0.nml[1]=0.0;
                    faces.push_back(fff0);
                    //      cnt++;                     
                }
            }                
        }
        //========
        //   cout<<cnt<<endl;
        current = current->next;
    }   
    cout<<facestmp.size()<<endl;
    facestmp.clear();
}


