#include"non_uniform_grid.h"

extern double dt;
extern Node *HeadListAllGrid;
extern OctCell *bodygrid;
extern double CFL;

void getMinTimeStep(OctCell *unp)
{	
	double ls;
	
	if( unp->flag == 0)
	{
		ls=unp->dt1;
		if(dt>ls) dt=ls; 
	}

}

void shblToJbbl(SHBL *, JBBL *);
void getTimeStep(OctCell *unp)
{

	double u, v;

	double dx,dy;

	double a;
	JBBL jbu;
	SHBL shu;

	dx=unp->dx;
	dy=unp->dy;
	shu.q=unp->dof[0][0];
	shu.qu=unp->dof[1][0];
	shu.qv=unp->dof[2][0];
	shu.te=unp->dof[3][0];	
			
	shblToJbbl(&shu,&jbu);

	u = jbu.u;
	v = jbu.v; 

	a = sqrt(GAMMA * jbu.p / jbu.q);

    	unp->dt1 = CFL * dx * dy/ ((fabs(u) + a)*dy + (fabs(v) + a)*dx);
}


void timestep()
{		
  Node *current;
  current = HeadListAllGrid;  
  OctCell *pcell0;
  
  dt=2000.0;			
  current = HeadListAllGrid;     
  while(current != NULL)
  {		
    pcell0 = current->cell;
    if(pcell0->flag == 0 && current->flg<=2 ) {
        getTimeStep(pcell0);
#ifndef LOCAL_TIME         
       getMinTimeStep(pcell0); 
#endif       
     }
    current = current->next;
  }			
}


