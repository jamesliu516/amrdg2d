
#include "non_uniform_grid.h"
#include"findneighbor.h"
//寻找物面相交的网格邻近的需要加密的网格单元地址
void findrefp(OctCell *pp[], OctCell *parent)
{                        		
/*	pp[0]=EastNeighbor(parent);
	pp[1]=NorthNeighbor(pp[0]);
       pp[2]=NorthNeighbor(parent);
	pp[3]=WestNeighbor(pp[2]);
	pp[4]=WestNeighbor(parent);
	pp[5]=SouthNeighbor(pp[4]);
	pp[6]=SouthNeighbor(parent);
	pp[7]=EastNeighbor(pp[6]);
	
	pp[8]=NorthNeighbor(pp[1]);

	pp[9]=WestNeighbor(pp[8]);
	pp[10]=WestNeighbor(pp[9]);
    pp[11]=WestNeighbor(pp[10]);
	pp[12]=WestNeighbor(pp[3]);
	pp[13]=WestNeighbor(pp[4]);
	pp[14]=WestNeighbor(pp[5]);
	pp[15]=SouthNeighbor(pp[14]);
	pp[16]=EastNeighbor(pp[15]);
	pp[17]=EastNeighbor(pp[16]);
	pp[18]=EastNeighbor(pp[17]);
	pp[19]=EastNeighbor(pp[18]);
    pp[20]=NorthNeighbor(pp[19]);
	pp[21]=NorthNeighbor(pp[20]);
	pp[22]=NorthNeighbor(pp[21]);
	pp[23]=NorthNeighbor(pp[22]);  */

	pp[0]=WestNeighbor(parent);
	pp[1]=WestNeighbor(pp[0]);
	pp[2]=WestNeighbor(pp[1]);
	pp[3]=EastNeighbor(parent);
	pp[4]=EastNeighbor(pp[3]);
	pp[5]=EastNeighbor(pp[4]);

	pp[6]=NorthNeighbor(pp[2]);
	pp[7]=NorthNeighbor(pp[6]);
	pp[8]=NorthNeighbor(pp[7]);

	pp[9]=SouthNeighbor(pp[2]);
	pp[10]=SouthNeighbor(pp[9]);
	pp[11]=SouthNeighbor(pp[10]);

	pp[12]=NorthNeighbor(pp[1]);
	pp[13]=NorthNeighbor(pp[12]);
	pp[14]=NorthNeighbor(pp[13]);

	pp[15]=SouthNeighbor(pp[1]);
	pp[16]=SouthNeighbor(pp[15]);
	pp[17]=SouthNeighbor(pp[16]);

	pp[18]=NorthNeighbor(pp[0]);
	pp[19]=NorthNeighbor(pp[18]);
	pp[20]=NorthNeighbor(pp[19]);

	pp[21]=SouthNeighbor(pp[0]);
	pp[22]=SouthNeighbor(pp[21]);
	pp[23]=SouthNeighbor(pp[22]);

	pp[24]=NorthNeighbor(parent);
	pp[25]=NorthNeighbor(pp[24]);
	pp[26]=NorthNeighbor(pp[25]);

	pp[27]=SouthNeighbor(parent);
	pp[28]=SouthNeighbor(pp[27]);
	pp[29]=SouthNeighbor(pp[28]);

	pp[30]=NorthNeighbor(pp[3]);
	pp[31]=NorthNeighbor(pp[30]);
	pp[32]=NorthNeighbor(pp[31]);

	pp[33]=SouthNeighbor(pp[3]);
	pp[34]=SouthNeighbor(pp[33]);
	pp[35]=SouthNeighbor(pp[34]);

	pp[36]=NorthNeighbor(pp[4]);
	pp[37]=NorthNeighbor(pp[36]);
	pp[38]=NorthNeighbor(pp[37]);

	pp[39]=SouthNeighbor(pp[4]);
	pp[40]=SouthNeighbor(pp[39]);
	pp[41]=SouthNeighbor(pp[40]);

	pp[42]=NorthNeighbor(pp[5]);
	pp[43]=NorthNeighbor(pp[42]);
	pp[44]=NorthNeighbor(pp[43]);

	pp[45]=SouthNeighbor(pp[5]);
	pp[46]=SouthNeighbor(pp[45]);
	pp[47]=SouthNeighbor(pp[46]);
}
