
#include"non_uniform_grid.h"
#include<ctime>

void out_para( const string &str123);
void outputcell(int ij=0);
void outputcell_new(int ij=0);
void outputcell_sol(int ij=0, const string& str123="sol");
void background_grid();
void pointsAtSolid();
void creatgrid();
void formListForAllGrid();
void setFlagForAll();
void face_vec();
void face_vec1();
void point_vec();
void formListFaceCellAll();
void releaseList();
void releaseFaceCellList();
void initialDataGrid();
void timestep();
void runge_kutta(int);
void refineCicleDemain(int rfN); // 1,2,3,4,5次加密
void SolutionAMR(int);
void initiallevel();
void find_multivalue_cell(const PXYZ &sharppoint);
extern OctCell *bodygrid;
extern double timesum;
extern double dt;
extern double CFL;

extern PXYZ MultiCircle[N_BODY+1][500]; //points at solid wall
//extern PXYZ Circle[NPOINTCIRCLE+1];

//const PXYZ *const wallpoint = Circle;
//const int NumberAirfoilPoint=NPOINTCIRCLE;
//const int HalfNumberAirfoilPoint=NumberAirfoilPoint / 2;

extern PXYZ ext_wall[NPT_EXT_WALL+1]; //points at external computational domain
void computeNormalArea();
void wherePointAll();
void shapeTransfer();
extern PXYZ MBsharp_point[][8];
void cal_te_res();
void output4restart();
void read4restart();
void form_inter_cell_set();
void formExtWallInOutBndry();
void set_exFlag4All(int nmpnt, PXYZ wallpt[]);//网格中心在固壁内部的为-1流场内部为0 在翼型上的为2, 尖后缘为4
extern double value_kxrcf_amr;
int main()
{
    clock_t start, finish; 
    double duration; 

    bool is_stop=false; 
    bool b_in_set=false;

    int nstp=0;
    int lscc=0; 
    int c_amr=1;
    string str123("_");
    double tpr0,tpr1;
    bool btpr0=true, btpr1=true;
    bool bl00=false,bl11=false;
    tpr0=0.08;
    tpr1=0.125;

    if(!(nDOF==3||nDOF==6)) {
        cout<<"We only consider p1 or p2 polynomial.\nSo nDOF=3 or nDOF=6.\nOthers are wrong.\n";
        cout<<"You input nDOF="<<nDOF<<".\n";
        cout<<"The code stop here."<<endl;
        return 0;
    }
    // cout<<"please input a string to add the file name: ";
    // cin>>str123;

    out_para(str123);
    bodygrid=new OctCell[Nx*Ny+1];

    background_grid();

    pointsAtSolid(); 
    shapeTransfer();
    creatgrid(); 

    refineCicleDemain(localRefine); // 1,2,3,4,5次加密  

    formListForAllGrid();  

    setFlagForAll();
    formListFaceCellAll();   
#if IS_TUBE==1
    set_exFlag4All(NPT_EXT_WALL, ext_wall);
    formExtWallInOutBndry();
#endif
    //  start = clock(); 
    face_vec();
    //   finish = clock();   
    //   duration = static_cast<double>(finish - start) / CLOCKS_PER_SEC;   
    //   cout<<duration<<" seconds-faces."<<endl; 

    initialDataGrid();
//    read4restart();
    initiallevel();
    //  start = clock();     
    point_vec();     
    outputcell_new(0);//output mesh
//    int j8;
  //  cin >>j8;
    //    finish = clock();  
    //   duration = static_cast<double>(finish - start) / CLOCKS_PER_SEC;   
    //   cout<<duration<<" seconds-output mesh modify.\n"<<endl;     

    //   start = clock();         
    //    outputcell(10);
    //    finish = clock();  
    //   duration = static_cast<double>(finish - start) / CLOCKS_PER_SEC;   
    //   cout<<duration<<" seconds-output mesh."<<endl;   
    
#if HAVE_MV!=0
    for(int wy=0; wy<N_BODY;++wy){
        for(int i_t=0;i_t<MVptNum[wy];++i_t){  
            MBsharp_point[wy][i_t]
                =MultiCircle[wy][n_sharp_point[wy][i_t]];

            find_multivalue_cell(MBsharp_point[wy][i_t]);
        }
    }
#endif
          
    computeNormalArea();
    form_inter_cell_set();
    wherePointAll();
    outputcell_sol(0,str123);
    // outputcell_new(0);
    for(int n=0; ; n++) {
#ifndef FIXED_DT
        timestep();
#else
        dt=fixed_dt;
#endif

#ifndef LOCAL_TIME         
        if(timesum+dt > t_print){
            dt=t_print-timesum;
            is_stop=true;
        }

        if(btpr0 && timesum+dt>tpr0){
            dt=tpr0-timesum;
            btpr0=false;
            bl00=true;
            CFL=0.034;
        }

        if(btpr1 && timesum+dt>tpr1){
            dt=tpr1-timesum;
            btpr1=false;
            bl11=true;
        }
#endif    
        //--------------------------------------------------------------
        runge_kutta(n);  
        if(n%20==0){
            cout<<"n: "<<n<<", dt="<<dt<<", time="<<timesum<<
                "  ";
            cal_te_res();
        }
        //---------------------------------------------------------------
#ifndef LOCAL_TIME       
        timesum+=dt;
        if( bl00 ) {			
            if(b_in_set){
                form_inter_cell_set();
                wherePointAll();
                b_in_set=false;
            }
            outputcell_sol(-time(NULL),"_t0p11");
            bl00=false;
        }

        if( bl11 ) {			
            if(b_in_set){
                form_inter_cell_set();
                wherePointAll();
                b_in_set=false;
            }
            outputcell_sol(-time(NULL),"_t0p2");
            bl11=false;
        }

        if(timesum>0.15)value_kxrcf_amr=1.0;
#endif
        if( (n%201 == 0&&n>=4000 )||(n<4000 && n%201==0&&n>0) ) {
//        if( (n%301 == 0 && n>0 ) ) {			
            lscc++;
            if(b_in_set){
                form_inter_cell_set();
                wherePointAll();
                b_in_set=false;
            }
            outputcell_sol(lscc,str123);
        }


     //   if(n==12000) output4restart();
        if(n%200==0 && n>=100 && SwAMR==true){
         //   if(n%1000==0 && n>=800 && SwAMR==true&& c_amr<=6){
    //    if(((n%200==0 && n>300 && n<2000)||(n%300==0&&n>2100) )
      //  if(((n%300==0&&n>200) )
        //    && SwAMR==true){
            SolutionAMR(c_amr);
            b_in_set=true;
            ++c_amr;
        }

        if(n>=max_time_steps&&n%500==0){
            cout<<"input an integer: 0 stop, !=0 continue to compute"<<endl;
            cin>>nstp;
            if(nstp==0) 
                is_stop=true;
        }

        if(is_stop) {
#ifndef LOCAL_TIME          
            cout<<"End time: "<<t_print<<endl;
#endif          
            if(b_in_set){
                form_inter_cell_set();
                wherePointAll();
                b_in_set=false;
            }
            outputcell_sol(-time(NULL));
            break;
        }
    }


    releaseList(); 
    releaseFaceCellList();
    delete [] bodygrid;     
    return 0;
    }
