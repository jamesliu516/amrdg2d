#############################################################
# DMU PhD research code
#16/11/2010, developed by Liu Jianming, Xuzhou Normal Univ.
#Version 0 
############################################################

CC=g++

OPTS= -O3
#OPTS= -Wall -g
EXE=pos_amrdg2

SRC=background_grid.cpp   bodygridclass.cpp     vec3d.cpp\
      commons.cpp  setflag.cpp   vec4d.cpp point_vec.cpp\
      outputgrid.cpp  findneighbor.cpp face_vec.cpp\
     points_solid.cpp     vecnd.cpp findrefp.cpp\
     formlist.cpp  amr_grid_treat.cpp   vec2d.cpp  main.cpp\
     approx_sol.cpp face_flux_all.cpp hllc.cpp roe.cpp\
     spatial_discretization.cpp wall4amrdg2013.cpp\
     rk.cpp farfield_boundary.cpp time_step.cpp\
     initial_data.cpp tvbm_limiter4all.cpp\
     tree_avg.cpp sonCellAvg.cpp out_para.cpp\
	 solveAMR.cpp multi_value_cell.cpp cal_te_res.cpp\
	  riemannWallBound.cpp

#wall4_amr_dg.cpp

.SUFFIXES: .cpp .o

OBJ=	$(SRC:.cpp=.o)

.cpp.o:
	$(CC) $(OPTS) -c $<

all:	$(EXE)

$(EXE):	$(OBJ)
	$(CC) $(OPTS) -o $@ $(OBJ)

clean:	
	rm  $(OBJ) $(EXE) 
