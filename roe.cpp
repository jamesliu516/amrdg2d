
#include"non_uniform_grid.h"

//void RoeRiemannSolver(JBBL *Ul, JBBL *Ur, JBBL *Uface, double *a) 

void RoeRiemannSolver(JBBL *Ul, JBBL *Ur, JBBL *Uface, double *a) 
{
	double ql,qr;

	double htl, htr,faceht;
	double dls;

	ql=sqrt(Ul->q);
	qr=sqrt(Ur->q);

 
	Uface->q = ql*qr;

	Uface->u=(ql*Ul->u + qr*Ur->u)/(ql+qr);
	Uface->v=(ql*Ul->v + qr*Ur->v)/(ql+qr);

        htl=Ul->p/(Ul->q * (GAMMA-1)) 
		+ 0.5 * (Ul->u * Ul->u + Ul->v * Ul->v)
		+ Ul->p/Ul->q;
        htr=Ur->p/(Ur->q * (GAMMA-1)) 
		+ 0.5 * (Ur->u * Ur->u + Ur->v * Ur->v)
		+ Ur->p/Ur->q;

	faceht=(ql*htl + qr*htr)/(ql+qr);

       htl=(faceht - 0.5 * (Uface->u*Uface->u +
		Uface->v*Uface->v)) * (GAMMA-1);

	Uface->p = Uface->q * htl  / GAMMA;
	*a=sqrt(htl);  // ÒôËÙ
}


