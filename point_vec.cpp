
#include"non_uniform_grid.h"
#include"findneighbor.h"
#include"vec2d.h"
#include<set>
extern Node *HeadListAllGrid;

extern OctCell *bodygrid;

extern map<OctCell *, cellPointIndex> cellPoints;      
extern vector<pointInCell> points4out; 	   

//相邻网格level差最多是1至关重要

void point_vec()
{
   Node *current;
   
   OctCell *pcell0;
   OctCell *pcell1;
   Vec2D x0,x1;
   
  // OctCell *pcellv[4]={NULL}; 
   map<OctCell *, short> pcellv; //共顶点的cell,第二个整数值表示cell的第几个顶点共此cell
   set<OctCell *> pcell_set;   //可能共点的那些cell
   
   cellPointIndex cpi1;	
   pointInCell pic1;
   
   int cnt=0;//the number of points
   short cnl=0;//表示点的类型有几个cell共此点,一般4还有自适应粗细界面处为2
   
   int idflg=0;
   
   int id1=0;
   short ii,jj;
   
   OctCell *pp[8]={NULL};
   OctCell *parent;
   
   current = HeadListAllGrid;  
   
   cellPoints.clear();
   points4out.clear();
  while(current != NULL)
  {		
    pcell0 = current->cell;
    if(pcell0->flag % 2 == 0 && current->flg==0) 
//    if(pcell0->exflag  == 2  ) 
         cellPoints.insert(make_pair(pcell0,cpi1)); 
    current = current->next;
  }
  
  cout<<"Total output cells: "<<cellPoints.size()<<endl;

  
  cnt=0;
  	   
    map<OctCell *, cellPointIndex>::iterator itmb;

    for(map<OctCell *, cellPointIndex>::iterator itmap=cellPoints.begin(); itmap!=cellPoints.end();++itmap) {
    
        parent=pcell0=itmap->first;       
    
        pp[0]=EastNeighbor(parent);
	 pp[1]=NorthNeighbor(pp[0]);
        pp[2]=NorthNeighbor(parent);
	 pp[3]=WestNeighbor(pp[2]);
	 pp[4]=WestNeighbor(parent);
	 pp[5]=SouthNeighbor(pp[4]);
	 pp[6]=SouthNeighbor(parent);
	 pp[7]=EastNeighbor(pp[6]);
         
        pcell_set.clear();
        if(pp[0]->reflag==0){
            pcell_set.insert(pp[0]);
         }         
        else {
            pcell_set.insert(pp[0]->children[3]);
            pcell_set.insert(pp[0]->children[0]);
        }
        
        if(pp[1]->reflag==0){
            pcell_set.insert(pp[1]);
         }         
        else {
            pcell_set.insert(pp[1]->children[0]);
         }        
         
        if(pp[2]->reflag==0){
            pcell_set.insert(pp[2]);
         }         
        else {
            pcell_set.insert(pp[2]->children[1]);
            pcell_set.insert(pp[2]->children[0]);
        }
        
        if(pp[3]->reflag==0){
            pcell_set.insert(pp[3]);
         }         
        else {
            pcell_set.insert(pp[3]->children[1]);
        }

        if(pp[4]->reflag==0){
            pcell_set.insert(pp[4]);
         }         
        else {
            pcell_set.insert(pp[4]->children[2]);
            pcell_set.insert(pp[4]->children[1]);
        }
                                
        if(pp[5]->reflag==0){
            pcell_set.insert(pp[5]);
         }         
        else {
            pcell_set.insert(pp[5]->children[2]);
        }
        
        if(pp[6]->reflag==0){
            pcell_set.insert(pp[6]);
         }         
        else {
            pcell_set.insert(pp[6]->children[3]);
            pcell_set.insert(pp[6]->children[2]);
        }
        
        if(pp[7]->reflag==0){
            pcell_set.insert(pp[7]);
         }         
        else {
            pcell_set.insert(pp[7]->children[3]);
        }        
        
         
        pcellv.clear(); 
        if(itmap->second.iPoint[0]==0) 
         {
             x0[0]=pcell0->xc1-0.5*pcell0->dx;
             x0[1]=pcell0->yc1-0.5*pcell0->dy;  
                                   
             idflg=0;  //取得已经标记的点对应点数组中的点序号  
             cnl=0;//表示点的类型有几个cell共此点,一般4还有自适应粗细界面处为2
             for(set<OctCell*>::iterator iter=pcell_set.begin(); iter!=pcell_set.end(); ++iter){
                  pcell1=*iter;
                  for(short i1=0;i1<4; ++i1){
 	              switch(i1) {	
	                  case 0: ii=-1;jj=-1;break;
	                  case 1: ii=1;jj=-1;break;
	                  case 2: ii=1;jj=1;break;
	                  case 3: ii=-1;jj=1;break;	
	                }            
	           x1[0]=pcell1->xc1+ii*pcell1->dx/2.0;				
	           x1[1]=pcell1->yc1+jj*pcell1->dy/2.0;	
	                
	           if(fabs(x1-x0)<ERRS){
	              ++cnl;
                     itmb=cellPoints.find(pcell1);                   
                     if(itmb!=cellPoints.end()){
                          pcellv.insert(make_pair(pcell1,i1));
                       }
                    }
                  }
                if(cnl==4) break;
             }
             
            for(map<OctCell*, short>::iterator itm1=pcellv.begin(); itm1!=pcellv.end(); ++itm1){
               if(cellPoints[itm1->first].iPoint[itm1->second]!=0)
                    idflg=cellPoints[itm1->first].iPoint[itm1->second];
             }
            
            if(idflg==0){
               cnt++;
               id1=cnt;
               pic1.pt=x0;  
               points4out.push_back(pic1);        
             }
            else{
               id1=idflg;
             }
             
           itmap->second.iPoint[0]=id1;
           
            for(map<OctCell*, short>::iterator itm1=pcellv.begin(); itm1!=pcellv.end(); ++itm1){
                    cellPoints[itm1->first].iPoint[itm1->second]=id1;
             }           
        } 
 //==          
        pcellv.clear(); 
        if(itmap->second.iPoint[1]==0) 
         {
             x0[0]=pcell0->xc1+0.5*pcell0->dx;
             x0[1]=pcell0->yc1-0.5*pcell0->dy;  
                                   
             idflg=0;  //取得已经标记的点对应点数组中的点序号  
             cnl=0;//表示点的类型有几个cell共此点,一般4还有自适应粗细界面处为3
             for(set<OctCell*>::iterator iter=pcell_set.begin(); iter!=pcell_set.end(); ++iter){
                  pcell1=*iter;
                  for(short i1=0;i1<4; ++i1){
 	              switch(i1) {	
	                  case 0: ii=-1;jj=-1;break;
	                  case 1: ii=1;jj=-1;break;
	                  case 2: ii=1;jj=1;break;
	                  case 3: ii=-1;jj=1;break;	
	                }            
	           x1[0]=pcell1->xc1+ii*pcell1->dx/2.0;				
	           x1[1]=pcell1->yc1+jj*pcell1->dy/2.0;	
	                
	           if(fabs(x1-x0)<ERRS){
	              ++cnl;
                     itmb=cellPoints.find(pcell1);                   
                     if(itmb!=cellPoints.end()){
                          pcellv.insert(make_pair(pcell1,i1));
                       }
                    }
                  }
                if(cnl==4) break;
             }
             
            for(map<OctCell*, short>::iterator itm1=pcellv.begin(); itm1!=pcellv.end(); ++itm1){
               if(cellPoints[itm1->first].iPoint[itm1->second]!=0)
                    idflg=cellPoints[itm1->first].iPoint[itm1->second];
             }
            
            if(idflg==0){
               cnt++;
               id1=cnt;
               pic1.pt=x0;                
                points4out.push_back(pic1);                        
             }
            else{
               id1=idflg;
             }
             
           itmap->second.iPoint[1]=id1;
           
            for(map<OctCell*, short>::iterator itm1=pcellv.begin(); itm1!=pcellv.end(); ++itm1){
                    cellPoints[itm1->first].iPoint[itm1->second]=id1;
             }           
        } 
//==
        pcellv.clear(); 
        if(itmap->second.iPoint[2]==0) 
         {
             x0[0]=pcell0->xc1+0.5*pcell0->dx;
             x0[1]=pcell0->yc1+0.5*pcell0->dy;  
                                   
             idflg=0;  //取得已经标记的点对应点数组中的点序号  
             cnl=0;//表示点的类型有几个cell共此点,一般4还有自适应粗细界面处为3
             for(set<OctCell*>::iterator iter=pcell_set.begin(); iter!=pcell_set.end(); ++iter){
                  pcell1=*iter;
                  for(short i1=0;i1<4; ++i1){
 	              switch(i1) {	
	                  case 0: ii=-1;jj=-1;break;
	                  case 1: ii=1;jj=-1;break;
	                  case 2: ii=1;jj=1;break;
	                  case 3: ii=-1;jj=1;break;	
	                }            
	           x1[0]=pcell1->xc1+ii*pcell1->dx/2.0;				
	           x1[1]=pcell1->yc1+jj*pcell1->dy/2.0;	
	                
	           if(fabs(x1-x0)<ERRS){
	              ++cnl;
                     itmb=cellPoints.find(pcell1);                   
                     if(itmb!=cellPoints.end()){
                          pcellv.insert(make_pair(pcell1,i1));
                       }
                    }
                  }
                if(cnl==4) break;
             }
             
            for(map<OctCell*, short>::iterator itm1=pcellv.begin(); itm1!=pcellv.end(); ++itm1){
               if(cellPoints[itm1->first].iPoint[itm1->second]!=0)
                    idflg=cellPoints[itm1->first].iPoint[itm1->second];
             }
            
            if(idflg==0){
               cnt++;
               id1=cnt;
               pic1.pt=x0;
               points4out.push_back(pic1);                
                         
             }
            else{
               id1=idflg;
             }
             
           itmap->second.iPoint[2]=id1;
           
            for(map<OctCell*, short>::iterator itm1=pcellv.begin(); itm1!=pcellv.end(); ++itm1){
                    cellPoints[itm1->first].iPoint[itm1->second]=id1;
             }           
        }         
// 

        pcellv.clear(); 
        if(itmap->second.iPoint[3]==0) 
         {
             x0[0]=pcell0->xc1-0.5*pcell0->dx;
             x0[1]=pcell0->yc1+0.5*pcell0->dy;  
                                   
             idflg=0;  //取得已经标记的点对应点数组中的点序号  
             cnl=0;//表示点的类型有几个cell共此点,一般4还有自适应粗细界面处为3
             for(set<OctCell*>::iterator iter=pcell_set.begin(); iter!=pcell_set.end(); ++iter){
                  pcell1=*iter;
                  for(short i1=0;i1<4; ++i1){
 	              switch(i1) {	
	                  case 0: ii=-1;jj=-1;break;
	                  case 1: ii=1;jj=-1;break;
	                  case 2: ii=1;jj=1;break;
	                  case 3: ii=-1;jj=1;break;	
	                }            
	           x1[0]=pcell1->xc1+ii*pcell1->dx/2.0;				
	           x1[1]=pcell1->yc1+jj*pcell1->dy/2.0;	
	                
	           if(fabs(x1-x0)<ERRS){
	              ++cnl;
                     itmb=cellPoints.find(pcell1);                   
                     if(itmb!=cellPoints.end()){
                          pcellv.insert(make_pair(pcell1,i1));
                       }
                    }
                  }
                if(cnl==4) break;
             }
             
            for(map<OctCell*, short>::iterator itm1=pcellv.begin(); itm1!=pcellv.end(); ++itm1){
               if(cellPoints[itm1->first].iPoint[itm1->second]!=0)
                    idflg=cellPoints[itm1->first].iPoint[itm1->second];
             }
            
            if(idflg==0){
               cnt++;
               id1=cnt;
               pic1.pt=x0; 
                points4out.push_back(pic1);                        
             }
            else{
               id1=idflg;
             }
             
           itmap->second.iPoint[3]=id1;
           
            for(map<OctCell*, short>::iterator itm1=pcellv.begin(); itm1!=pcellv.end(); ++itm1){
                    cellPoints[itm1->first].iPoint[itm1->second]=id1;
             }           
        }        
   }     
   
   cout<<"Total Output Points: "<<cnt<<endl; 
   bool is_here=false;
   int tnmc;
   for(map<OctCell *, cellPointIndex>::iterator itmap=cellPoints.begin(); itmap!=cellPoints.end();++itmap) {
        pcell0=itmap->first;
       for(int m=0; m<4;++m) {
         is_here=false;
         tnmc=points4out[itmap->second.iPoint[m]-1].nmc;
         for(int i1=0;i1<tnmc; ++i1) {
             if(points4out[itmap->second.iPoint[m]-1].mcell[i1]==pcell0) {
                  is_here=true;
                  break;
               }
          }
         if(!is_here) {
              points4out[itmap->second.iPoint[m]-1].mcell[tnmc]=pcell0;
              points4out[itmap->second.iPoint[m]-1].nmc+=1;
         }
      }
    }
    
  }                
                

