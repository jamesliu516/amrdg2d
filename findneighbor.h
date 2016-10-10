#ifndef FIND_NEIGHBOR_H
#define FIND_NEIGHBOR_H
#include "non_uniform_grid.h"

OctCell *NorthNeighbor(OctCell *pp);
OctCell *WestNeighbor(OctCell *pp);
OctCell *EastNeighbor(OctCell *pp);
OctCell *SouthNeighbor(OctCell *pp);

#endif
