
#include"non_uniform_grid.h"

//extern PXYZ Circle[NPOINTCIRCLE+1];

extern PXYZ MultiCircle[N_BODY+1][500]; //points at solid wall
extern int NWallPts[N_BODY+3];
extern string WallPtsFile[N_BODY+3];

//const int NumberAirfoilPoint=NPOINTCIRCLE;

//const PXYZ *const wallpoint = Circle;
//const int HalfNumberAirfoilPoint=NumberAirfoilPoint / 2;

extern PXYZ ext_wall[NPT_EXT_WALL+1]; //points at external computational domain
void pointsAtSolid()
{
    int i;

    string filename;
    double stp; 
    double x00,y00;
    ifstream infile;
    for(int jj=0;jj<N_BODY;++jj) {

        filename=WallPtsFile[jj];

        cout << filename.c_str() << endl;
        infile.open(filename.c_str());

        stp = 2.0*PI/NWallPts[jj];
        if (jj==0){
            x00=0.0;
            y00=0.0;
        }
        else if(jj==1){
            x00=-1.5;
            y00=-0.50;
        }

        for(i=0; i < NWallPts[jj]; i++)
        {	
            infile >>MultiCircle[jj][i].x;
            infile >> MultiCircle[jj][i].y;
       //  if(jj==0)   MultiCircle[jj][i].x=MultiCircle[0][i].x+0.03;
//            if(jj==1) MultiCircle[jj][i].x=MultiCircle[0][i].x+0.5;
  //          if(jj==1) MultiCircle[jj][i].y=MultiCircle[0][i].y+0.5;
//              MultiCircle[jj][i].x =x00+ RADII * cos(PI+i*stp);
  //           MultiCircle[jj][i].y =y00+RADII * sin(PI+i*stp);
            //	cout<<Circle[i].x<<"  "<<Circle[i].y<<endl;
        }

        MultiCircle[jj][NWallPts[jj]]=MultiCircle[jj][0];
        infile.close();
        infile.clear();
    }
#if IS_TUBE ==1
    filename=ExtWALL_POINT_FILE;

    cout << filename.c_str() << endl;
    infile.open(filename.c_str());

    for(i=0; i < NPT_EXT_WALL; i++)
    {	
        infile >> ext_wall[i].x;
        infile >> ext_wall[i].y;
    }
    ext_wall[NPT_EXT_WALL]=ext_wall[0];

    infile.close();
    infile.clear();
#endif
}

void shapeTransfer()
{
    int i;
    double phi, x1,y1, cphi,sphi;
    phi = rotate_deg * PI/180.0;
    cphi = cos(-phi);
    sphi = sin(-phi);

    for(int jj=0;jj<N_BODY;++jj) {
        for(i=0;i< NWallPts[jj]; i++)
        {
            x1 = MultiCircle[jj][i].x *  cphi
                - MultiCircle[jj][i].y * sphi;
            y1 = MultiCircle[jj][i].x *  sphi
                + MultiCircle[jj][i].y * cphi;
            MultiCircle[jj][i].x = x1;
            MultiCircle[jj][i].y = y1;
        }

        MultiCircle[jj][NWallPts[jj]]=MultiCircle[jj][0];
    }
}


int newFlagpoint(double x,double y, int NumberPoint,PXYZ wallpnt[] )
{
    int ie, cnt=0;
    PXYZ pnt1,pt1,pt2;
	pnt1.x=x;
	pnt1.y=y;
    for(ie=0; ie <= NumberPoint - 1;++ie){
        pt1.x=wallpnt[ie].x;
        pt1.y=wallpnt[ie].y;
        pt2.x=wallpnt[ie+1].x;
        pt2.y=wallpnt[ie+1].y;

        if ((pt1.y!=pt2.y)&&(((pnt1.y >= pt1.y) && (pnt1.y < pt2.y)) 
                || ((pnt1.y >= pt2.y) && (pnt1.y < pt1.y)))) 
        {
            if(pnt1.x == (pt2.x - pt1.x) * (pnt1.y - pt1.y) / (pt2.y - pt1.y) + pt1.x)
            {	
                cnt=-1;break;
            }
            else  if (pnt1.x < (pt2.x - pt1.x) * (pnt1.y - pt1.y) / (pt2.y - pt1.y) + pt1.x)
            {	
                cnt++;
            }
        }
    }          
    return (cnt%2>0)?-1:0;
}

int newFlagpoint_1(double x, double y, int NumberPoint,PXYZ wallpnt[] )//应用射线方法 //in the solid flag==-1, out the solid flag=0
    //at the surface of the solid flag=2
{
    int j;
    double x1, y1, x2, y2;
    double dls1, dls2, nx, ny, pax, pay , pbx, pby, px, py, kslop, bs;
    int flag, xDirflag=0,yDirflag=0;
    const double smallp = 1e-10;
    x1=x;
    y1=20.2;
    x2=20.2;
    y2=y;

    for(j=0; j <= NumberPoint - 1; j++)
    {
        if(j == NumberPoint - 1)
        {
            if(fabs(wallpnt[j].y - wallpnt[0].y) < smallp )
            {
                if((wallpnt[0].x - x) * (wallpnt[j].x - x) +1.0 <= 1.0)
                {
                    if(fabs(wallpnt[j].y - y) < smallp )
                    { 
                        flag=2; 
                        return flag;
                    }
                }
            }
            else if(fabs(wallpnt[j].x - wallpnt[0].x) < smallp )
            {
                if((wallpnt[0].y - y) * (wallpnt[j].y - y) +1.0<= 1.0)
                {
                    if(fabs(wallpnt[j].x - x) < smallp )
                    { 
                        flag=2; 
                        return flag;
                    }
                }
            }
            else
            {
                kslop = (wallpnt[0].y-wallpnt[j].y)/(wallpnt[0].x-wallpnt[j].x);
                bs = wallpnt[j].y - kslop * wallpnt[j].x; 
                if((wallpnt[0].x - x) * (wallpnt[j].x - x) +1.0<= 1.0 
                    && fabs(y - kslop * x - bs) < smallp ) 
                {
                    flag=2;
                    return flag;
                }
            }

            dls1 = (wallpnt[j].x - x1) * (wallpnt[0].x - x1);

            nx = wallpnt[0].y - wallpnt[j].y;
            ny = -(wallpnt[0].x - wallpnt[j].x);
            px = 0.5 * (wallpnt[0].x + wallpnt[j].x);
            py = 0.5 * (wallpnt[0].y + wallpnt[j].y);
        }
        else
        { 
            if(fabs(wallpnt[j].y - wallpnt[j+1].y) < smallp )
            {
                if((wallpnt[j+1].x - x) * (wallpnt[j].x - x) +1.0<= 1.0)
                {
                    if(fabs(wallpnt[j].y - y) < smallp )
                    { 
                        flag=2; 
                        return flag;
                    }
                }
            }
            else if(fabs(wallpnt[j].x - wallpnt[j+1].x) < smallp )
            {
                if((wallpnt[j+1].y - y) * (wallpnt[j].y - y) +1.0<= 1.0)
                {
                    if(fabs(wallpnt[j].x - x) < smallp )
                    { 
                        flag=2; 
                        return flag;
                    }
                }
            }
            else
            {
                kslop = (wallpnt[j+1].y-wallpnt[j].y)/(wallpnt[j+1].x-wallpnt[j].x);
                bs = wallpnt[j].y - kslop * wallpnt[j].x; 
                if((wallpnt[j+1].x - x) * (wallpnt[j].x - x) +1.0<= 1.0 
                    && fabs(y - kslop * x - bs) < smallp ) 
                {
                    flag=2;
                    return flag;
                }
            }

            dls1 = (wallpnt[j].x - x1) * (wallpnt[j+1].x - x1);

            nx = wallpnt[j+1].y - wallpnt[j].y;
            ny = -(wallpnt[j+1].x - wallpnt[j].x);
            px = 0.5 * (wallpnt[j+1].x + wallpnt[j].x);
            py = 0.5 * (wallpnt[j+1].y + wallpnt[j].y);
        }

        pax = x - px;
        pay = y - py;
        pbx = x1 - px;
        pby = y1 - py;

        dls2 = (nx * pax + ny * pay) * (nx * pbx + ny * pby);

        if (dls1 <= 0.0 && dls2 <= 0.0) yDirflag++;
    }

    if(yDirflag % 2 == 0) flag = 0;
    else flag = -1;
    return flag;
}
//+++++++++++++++++++++++++++++++==after here useless
/*
int SolveRayTracing(double x, double y,struct edge wallface[],int Wp)
{	           
	int ie,cnt=0;
	struct point pnt1,pt1,pt2;
	pnt1.x=x;
	pnt1.y=y;
	for(ie=0;ie<Wp;ie++)
	{
		pt1.x=wallface[ie].epoint[0]->x;pt1.y=wallface[ie].epoint[0]->y;
		pt2.x=wallface[ie].epoint[1]->x;pt2.y=wallface[ie].epoint[1]->y; 
		if ((pt1.y!=pt2.y)&&(((pnt1.y >= pt1.y) && (pnt1.y < pt2.y)) 
			|| ((pnt1.y >= pt2.y) && (pnt1.y < pt1.y)))) 
		{
			if(pnt1.x == (pt2.x - pt1.x) * (pnt1.y - pt1.y) / (pt2.y - pt1.y) + pt1.x)
			{	
				cnt=1;break;
			}
			else  if (pnt1.x < (pt2.x - pt1.x) * (pnt1.y - pt1.y) / (pt2.y - pt1.y) + pt1.x)
			{	
				cnt++;
			}
		}
	}          
	return (cnt%2>0)?1:0;
}
*/
/*
int newFlagpoint(double x,double y)
{
    int ie, cnt=0;
    PXYZ pnt1,pt1,pt2;
    pnt1.x=x;
    pnt1.y=y;
    for(ie=0; ie <= NumberAirfoilPoint - 1; ++ie){
        pt1.x=wallpoint[ie].x;
        pt1.y=wallpoint[ie].y;
        pt2.x=wallpoint[ie+1].x;
        pt2.y=wallpoint[ie+1].y;

        if ((pt1.y!=pt2.y)&&(((pnt1.y >= pt1.y) && (pnt1.y < pt2.y)) 
                || ((pnt1.y >= pt2.y) && (pnt1.y < pt1.y)))) 
        {
            if(pnt1.x == (pt2.x - pt1.x) * (pnt1.y - pt1.y) / (pt2.y - pt1.y) + pt1.x)
            {	
                cnt=-1;break;
            }
            else  if (pnt1.x < (pt2.x - pt1.x) * (pnt1.y - pt1.y) / (pt2.y - pt1.y) + pt1.x)
            {	
                cnt++;
            }
        }
    }          
    return (cnt%2>0)?-1:0;
}

int newFlagpoint_1(double x, double y)//应用射线方法 //in the solid flag==-1, out the solid flag=0
                                         //at the surface of the solid flag=2
{
	int j;
	double x1, y1, x2, y2;
	double dls1, dls2, nx, ny, pax, pay , pbx, pby, px, py, kslop, bs;
	int flag, xDirflag=0,yDirflag=0;
	const double smallp = 1e-10;
//	const PXYZ *wallpoint = Circle;
	x1=x;
	y1=20.2;
	x2=20.2;
	y2=y;

//	if(x < RXMIN - 1.0e-10 || x > RXMAX + 1.0e-10 || y<= RYMIN || y>=RYMAX) 
//	{
//		flag=0;
//		return flag;
//	}
//	else
	//{
		for(j=0; j <= NumberAirfoilPoint - 1; j++)
		{
			if(j == NumberAirfoilPoint - 1)
			{
				if(fabs(wallpoint[j].y - wallpoint[0].y) < smallp )
				{
					if((wallpoint[0].x - x) * (wallpoint[j].x - x) +1.0 <= 1.0)
					{
						if(fabs(wallpoint[j].y - y) < smallp )
						{ 
							flag=2; 
							return flag;
						}
					}
				}
				else if(fabs(wallpoint[j].x - wallpoint[0].x) < smallp )
				{
					if((wallpoint[0].y - y) * (wallpoint[j].y - y) +1.0<= 1.0)
					{
						if(fabs(wallpoint[j].x - x) < smallp )
						{ 
							flag=2; 
					     	return flag;
						}
					}
				}
				else
				{
					kslop = (wallpoint[0].y-wallpoint[j].y)/(wallpoint[0].x-wallpoint[j].x);
					bs = wallpoint[j].y - kslop * wallpoint[j].x; 
					if((wallpoint[0].x - x) * (wallpoint[j].x - x) +1.0<= 1.0 
						&& fabs(y - kslop * x - bs) < smallp ) 
					{
						flag=2;
						return flag;
					}
				}
		        
				dls1 = (wallpoint[j].x - x1) * (wallpoint[0].x - x1);

	     	   	nx = wallpoint[0].y - wallpoint[j].y;
		       	ny = -(wallpoint[0].x - wallpoint[j].x);
		        px = 0.5 * (wallpoint[0].x + wallpoint[j].x);
		       	py = 0.5 * (wallpoint[0].y + wallpoint[j].y);
			}
			else
			{ 
				if(fabs(wallpoint[j].y - wallpoint[j+1].y) < smallp )
				{
					if((wallpoint[j+1].x - x) * (wallpoint[j].x - x) +1.0<= 1.0)
					{
						if(fabs(wallpoint[j].y - y) < smallp )
						{ 
							flag=2; 
							return flag;
						}
					}
				}
				else if(fabs(wallpoint[j].x - wallpoint[j+1].x) < smallp )
				{
					if((wallpoint[j+1].y - y) * (wallpoint[j].y - y) +1.0<= 1.0)
					{
						if(fabs(wallpoint[j].x - x) < smallp )
						{ 
							flag=2; 
					     	return flag;
						}
					}
				}
				else
				{
					kslop = (wallpoint[j+1].y-wallpoint[j].y)/(wallpoint[j+1].x-wallpoint[j].x);
					bs = wallpoint[j].y - kslop * wallpoint[j].x; 
					if((wallpoint[j+1].x - x) * (wallpoint[j].x - x) +1.0<= 1.0 
						&& fabs(y - kslop * x - bs) < smallp ) 
					{
						flag=2;
						return flag;
					}
				}
		        
				dls1 = (wallpoint[j].x - x1) * (wallpoint[j+1].x - x1);

	     	   	nx = wallpoint[j+1].y - wallpoint[j].y;
		       	ny = -(wallpoint[j+1].x - wallpoint[j].x);
		        px = 0.5 * (wallpoint[j+1].x + wallpoint[j].x);
		       	py = 0.5 * (wallpoint[j+1].y + wallpoint[j].y);
			}

			pax = x - px;
			pay = y - py;
			pbx = x1 - px;
			pby = y1 - py;

			dls2 = (nx * pax + ny * pay) * (nx * pbx + ny * pby);

			if (dls1 <= 0.0 && dls2 <= 0.0) yDirflag++;
		}

		if(yDirflag % 2 == 0) flag = 0;
		else flag = -1;
		return flag;
	//}
}

*/
