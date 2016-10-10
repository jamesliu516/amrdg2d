
#include"non_uniform_grid.h"
#include"vec4d.h"

//void shblToJbbl(SHBL *, JBBL *);
//void jbblToShbl(JBBL *, SHBL *);
void RoeRiemannSolver(JBBL *Ul, JBBL *Ur, JBBL *Uface, double *a);

//计算通量函数的子程序
//HLLC Flux j.c.p 137, 38-78 (1997), or siam j.sci.comput. vol18, no6,1553-1570,1997
//void HLLC_flux(JBBL *ujbL, JBBL *ujbR, PXYZ *nml, Flux *flxface)
void HLLC_flux(JBBL *ujbL, JBBL *ujbR, PXYZ *nml, double flxface[4])
//sface 面积向量
{
	double s_l, s_r, s_m;
	double pStar;
	double aFace, ql, qr, al, ar, qFace; // speed

	double omgL, omgR;
//

	JBBL ul, ur, uFace;
	SHBL wLstar, wRstar, wl, wr;

	PXYZ normal;

	//Flux flxL, flxR, flxLs, flxRs;
	
	normal.x = nml->x;
	normal.y = nml->y;


	ul.q = ujbL->q;
	ul.u = ujbL->u;
	ul.v = ujbL->v;
	ul.p = ujbL->p;

	jbblToShbl(&ul, &wl);

	al = sqrt(GAMMA*ul.p/ul.q);

	ql = normal.x * ul.u + normal.y * ul.v;

	ur.q = ujbR->q;
	ur.u = ujbR->u;
	ur.v = ujbR->v;
	ur.p = ujbR->p;


	jbblToShbl(&ur, &wr);

	ar = sqrt(GAMMA*ur.p/ur.q);
	qr = normal.x * ur.u + normal.y * ur.v;
//
//
	RoeRiemannSolver(&ul, &ur, &uFace, &aFace);
	qFace = normal.x * uFace.u + normal.y * uFace.v;


	s_l = min2(ql-al, qFace-aFace);
	s_r = max2(qr+ar, qFace+aFace);


	s_m = (ur.q * qr * (s_r - qr) - ul.q * ql * (s_l - ql) + ul.p - ur.p)
			/ (ur.q * (s_r - qr) - ul.q * (s_l - ql));

	omgL = 1.0 / (s_l-s_m);
	omgR = 1.0 / (s_r-s_m);

	pStar = ul.q * (ql - s_l) * (ql - s_m) + ul.p;

	wLstar.q = omgL * ul.q * (s_l - ql);
	wLstar.qu = omgL * ((s_l - ql) * wl.qu + (pStar - ul.p) * normal.x);
	wLstar.qv = omgL * ((s_l - ql) * wl.qv + (pStar - ul.p) * normal.y);
	wLstar.te = omgL * ((s_l - ql) * wl.te - ul.p * ql + pStar * s_m);

	wRstar.q = omgR * ur.q * (s_r - qr);
	wRstar.qu = omgR * ((s_r - qr) * wr.qu + (pStar - ur.p) * normal.x);
	wRstar.qv = omgR * ((s_r - qr) * wr.qv + (pStar - ur.p) * normal.y);
	wRstar.te = omgR * ((s_r - qr) * wr.te - ur.p * qr + pStar * s_m);

	if(s_l > 0.0)
	{
		flxface[0] = ul.q * ql;
		flxface[1] = wl.qu * ql + ul.p * normal.x;
		flxface[2] = wl.qv * ql + ul.p * normal.y;
		flxface[3] = (wl.te + ul.p) * ql;


	//	flxface->fq = flxL.fq;
	//	flxface->fqu = flxL.fqu;
	//	flxface->fqv = flxL.fqv;
	//	flxface->fte = flxL.fte;

	}
	else if(s_r < 0.0)
	{
		flxface[0] = ur.q * qr;
		flxface[1] = wr.qu * qr + ur.p * normal.x;
		flxface[2] = wr.qv * qr + ur.p * normal.y;
		flxface[3] = (wr.te + ur.p) * qr;


	//	flxface->fq = flxR.fq;
	//	flxface->fqu = flxR.fqu;
	//	flxface->fqv = flxR.fqv;
	//	flxface->fte = flxR.fte;
	}
	else if(s_l <= 0.0 && s_m > 0.0)
	{
		flxface[0] = wLstar.q * s_m;
		flxface[1] = wLstar.qu * s_m + pStar * normal.x;
		flxface[2] = wLstar.qv * s_m + pStar * normal.y;
		flxface[3] = (wLstar.te + pStar) * s_m;

	//	flxface->fq = flxLs.fq;
	//	flxface->fqu = flxLs.fqu;
	//	flxface->fqv = flxLs.fqv;
	//	flxface->fte = flxLs.fte;

	}
	else if(s_m <= 0.0 && s_r >= 0.0)
	{
		flxface[0] = wRstar.q * s_m;
		flxface[1] = wRstar.qu * s_m + pStar * normal.x;
		flxface[2] = wRstar.qv * s_m + pStar * normal.y;
		flxface[3] = (wRstar.te + pStar) * s_m;

		//flxface->fq = flxRs.fq;
		//flxface->fqu = flxRs.fqu;
		//flxface->fqv = flxRs.fqv;
		//flxface->fte = flxRs.fte;	

	}
}

void HLLC_flux(SHBL &wl, SHBL &wr, PXYZ &normal, double flxface[4])
//sface 面积向量
{
	double s_l, s_r, s_m;
	double pStar;
	double aFace, ql, qr, al, ar, qFace; // speed

	double omgL, omgR;
//
	JBBL ul, ur, uFace;
	SHBL wLstar, wRstar;
	
	shblToJbbl(&wl, &ul);

	al = sqrt(GAMMA*ul.p/ul.q);

	ql = normal.x * ul.u + normal.y * ul.v;

	shblToJbbl(&wr, &ur);

	ar = sqrt(GAMMA*ur.p/ur.q);
	qr = normal.x * ur.u + normal.y * ur.v;
//
//
	RoeRiemannSolver(&ul, &ur, &uFace, &aFace);
	qFace = normal.x * uFace.u + normal.y * uFace.v;


	s_l = min2(ql-al, qFace-aFace);
	s_r = max2(qr+ar, qFace+aFace);


	s_m = (ur.q * qr * (s_r - qr) - ul.q * ql * (s_l - ql) + ul.p - ur.p)
			/ (ur.q * (s_r - qr) - ul.q * (s_l - ql));

	omgL = 1.0 / (s_l-s_m);
	omgR = 1.0 / (s_r-s_m);

	pStar = ul.q * (ql - s_l) * (ql - s_m) + ul.p;

	wLstar.q = omgL * ul.q * (s_l - ql);
	wLstar.qu = omgL * ((s_l - ql) * wl.qu + (pStar - ul.p) * normal.x);
	wLstar.qv = omgL * ((s_l - ql) * wl.qv + (pStar - ul.p) * normal.y);
	wLstar.te = omgL * ((s_l - ql) * wl.te - ul.p * ql + pStar * s_m);

	wRstar.q = omgR * ur.q * (s_r - qr);
	wRstar.qu = omgR * ((s_r - qr) * wr.qu + (pStar - ur.p) * normal.x);
	wRstar.qv = omgR * ((s_r - qr) * wr.qv + (pStar - ur.p) * normal.y);
	wRstar.te = omgR * ((s_r - qr) * wr.te - ur.p * qr + pStar * s_m);

	if(s_l > 0.0)
	{
		flxface[0] = ul.q * ql;
		flxface[1] = wl.qu * ql + ul.p * normal.x;
		flxface[2] = wl.qv * ql + ul.p * normal.y;
		flxface[3] = (wl.te + ul.p) * ql;

	}
	else if(s_r < 0.0)
	{
		flxface[0] = ur.q * qr;
		flxface[1] = wr.qu * qr + ur.p * normal.x;
		flxface[2] = wr.qv * qr + ur.p * normal.y;
		flxface[3] = (wr.te + ur.p) * qr;
	}
	else if(s_l <= 0.0 && s_m > 0.0)
	{
		flxface[0] = wLstar.q * s_m;
		flxface[1] = wLstar.qu * s_m + pStar * normal.x;
		flxface[2] = wLstar.qv * s_m + pStar * normal.y;
		flxface[3] = (wLstar.te + pStar) * s_m;
	}
	else if(s_m <= 0.0 && s_r >= 0.0)
	{
		flxface[0] = wRstar.q * s_m;
		flxface[1] = wRstar.qu * s_m + pStar * normal.x;
		flxface[2] = wRstar.qv * s_m + pStar * normal.y;
		flxface[3] = (wRstar.te + pStar) * s_m;
	}
}

extern double max_cv_x;
extern double max_cv_y;

void GLF_flux(SHBL &wl, SHBL &wr, PXYZ &nml, double flxface[4])  
{
    JBBL ul,ur;
    double vnl,vnr;
    double tzmax;
    
	shblToJbbl(&wl, &ul);
	shblToJbbl(&wr, &ur);

    vnl=ul.u*nml.x+ul.v*nml.y;
    vnr=ur.u*nml.x+ur.v*nml.y;
    if(fabs(nml.x)<1e-7) tzmax=max_cv_y;
    if(fabs(nml.y)<1e-7) tzmax=max_cv_x;

    flxface[0]=0.5*(wl.q*vnl+wr.q*vnr-tzmax*(wr.q-wl.q));
    flxface[1]=0.5*(wl.qu*vnl+nml.x*ul.p
        +wr.qu*vnr+nml.x*ur.p-tzmax*(wr.qu-wl.qu));
    flxface[2]=0.5*(wl.qv*vnl+nml.y*ul.p
        +wr.qv*vnr+nml.y*ur.p-tzmax*(wr.qv-wl.qv));

    flxface[3]=0.5*((wl.te+ul.p)*vnl + (wr.te+ur.p)*vnr
        - tzmax*(wr.te-wl.te));
}

extern Node *HeadListAllGrid;

void set_max_cv()
{
    Node *current;
    current = HeadListAllGrid;  
    OctCell *unp;
    double a;
    JBBL jbu;
    SHBL shu;
    double rt1,rt0;

    max_cv_x=max_cv_y=0.0;

    while(current != NULL)
    {		
        unp = current->cell;
        if(unp->flag  == 0 && current->flg<=2 ) {

            shu.q=unp->dof[0][0];
            shu.qu=unp->dof[1][0];
            shu.qv=unp->dof[2][0];
            shu.te=unp->dof[3][0];	

            shblToJbbl(&shu,&jbu);

            a = sqrt(GAMMA * jbu.p / jbu.q);
            rt0=fabs(jbu.u)+a;
            rt1=fabs(jbu.v)+a;
            if(max_cv_x<rt0)max_cv_x=rt0;
            if(max_cv_y<rt1)max_cv_y=rt1;
        }
        current = current->next;
    }
}




