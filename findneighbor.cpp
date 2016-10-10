#include "non_uniform_grid.h"

extern OctCell *bodygrid;
//ураз╬с
OctCell *NorthNeighbor(OctCell *pp)
{
	int NOchildren;
	long NOparent;
	OctCell *parent,*p;
    if(pp==NULL){
        return NULL;
    }
    else{
	parent=pp->parent;
	NOparent=pp->NOparent;
	if(parent==NULL) 
	{
		return(&bodygrid[NOparent+Nx]);
	}
	
	NOchildren=pp->NOchildren;
	if(NOchildren==0) return(parent->children[3]);
    if(NOchildren==1) return(parent->children[2]);
    p=NorthNeighbor(parent);
	if(p->children[0]==NULL) return(p);
	else if(NOchildren==3) return(p->children[0]);
	else if(NOchildren==2) return(p->children[1]);
    }
}

OctCell *WestNeighbor(OctCell *pp)
{
	int NOchildren;
	long NOparent;
	OctCell *parent,*p;
    if(pp==NULL){
        return NULL;
    }
    else{
	parent=pp->parent;
	NOparent=pp->NOparent;

	if(parent==NULL) return(&bodygrid[NOparent-1]);	

	NOchildren=pp->NOchildren;
	if(NOchildren==1) return(parent->children[0]);
    if(NOchildren==2) return(parent->children[3]);
    p=WestNeighbor(parent);
	if(p->children[0]==NULL) return(p);
	else if(NOchildren==0) return(p->children[1]);
	else if(NOchildren==3) return(p->children[2]);
    }
}

OctCell *EastNeighbor(OctCell *pp)
{
	int NOchildren;
	long NOparent;
	OctCell *parent,*p;
    if(pp==NULL){
        return NULL;
    }
    else{
	parent=pp->parent;
	NOparent=pp->NOparent;
	if(parent==NULL) return(&bodygrid[NOparent+1]);	
	NOchildren=pp->NOchildren;
	if(NOchildren==0) return(parent->children[1]);
    if(NOchildren==3) return(parent->children[2]);
    p=EastNeighbor(parent);
	if(p->children[0]==NULL) return(p);
	else if(NOchildren==1) return(p->children[0]);
	else if(NOchildren==2) return(p->children[3]);
    }
}

OctCell *SouthNeighbor(OctCell *pp)
{
	int NOchildren;
	long NOparent;
	OctCell *parent,*p;
    if(pp==NULL){
        return NULL;
    }
    else{
	parent=pp->parent;
	NOparent=pp->NOparent;
	if(parent==NULL) return(&bodygrid[NOparent-Nx]);	
	NOchildren=pp->NOchildren;
	if(NOchildren==3) return(parent->children[0]);
    if(NOchildren==2) return(parent->children[1]);
    p=SouthNeighbor(parent);
	if(p->children[0]==NULL) return(p);
	else if(NOchildren==0) return(p->children[3]);
	else if(NOchildren==1) return(p->children[2]);
    }
}
