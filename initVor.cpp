
void initialvortex(struct cell *pc,double U[])
{
	double dtT,dU,tau,rxy,epsi=0.5,Radc=0.05,alph=0.204;
	double x,y,xo,yo;
	xo=0.25,yo=0.5;
	x=pc->xc;y=pc->yc;
//	if(x<0.5)
	{					
		rxy=(x-xo)*(x-xo)+(y-yo)*(y-yo);
		if(rxy>1e-6)rxy=sqrt(rxy);

		tau=rxy/Radc;
		dtT=-(gama-1.0)*epsi*epsi*exp(2.0*alph*(1.0-tau*tau))/(4.0*alph*gama);
		dU=epsi*tau*exp(alph*(1.0-tau*tau));

		U[0]=pow(1.0+dtT,2.5);
		if(rxy>1e-6){
			U[1]=Uinf[1]+dU*(y-yo)/rxy;
			U[2]=Uinf[2]-dU*(x-xo)/rxy;
		}else{		
			U[1]=Uinf[1];	
			U[2]=Uinf[2];
		}
		U[3]=Uinf[3];
		U[4]=pow(U[0],1.4);			
	}
/*	else				
	{					
		U[0]=1.169;				
		U[1]=1.1133;					
		U[2]=0.0;
		U[3]=0.0;
		U[4]=1.245;				
	}*/
//	printf("d=%f,u=%f,v=%f,w=%f,p=%f\n",U[0],U[1],U[2],U[3],U[4]);
}
