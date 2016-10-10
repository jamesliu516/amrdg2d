#include <fstream>
#include <vector>
#include <cmath>
#include<numeric>
using namespace std;

#include "vecnd.h"

//Vec3D::Vec3D(double x, double y, double z) {
//	comp[0]=x;
//	comp[1]=y;
//	comp[2]=z;
//}

VecND::VecND(int n) {
    comp.resize(n,0.0);
}

void VecND::resize(int n, double r) {
         comp.clear();
         comp.resize(n,r);
}


int VecND::size() {
    return comp.size();
}

double VecND::dot(const VecND &right) {
    //   double rtmp=0.0;
    //   for(vdsz_type cnt=0; cnt!=comp.size(); ++cnt)
   //           rtmp+=comp[cnt]*right.comp[cnt];
   //    return rtmp;
    return (inner_product(comp.begin(), comp.end(), right.comp.begin(),0.0));
//	return (comp[0]*right.comp[0]+comp[1]*right.comp[1]+comp[2]*right.comp[2]);
}

double fabs(const VecND &vec) {
	//return sqrt(vec.comp[0]*vec.comp[0]+vec.comp[1]*vec.comp[1]+vec.comp[2]*vec.comp[2]);
	//return sqrt(vec.dot(vec));
	
       double rtmp=0.0;
       for(vdsz_type cnt=0; cnt!=vec.comp.size(); ++cnt)
              rtmp+=vec.comp[cnt]*vec.comp[cnt];
       return sqrt(rtmp);	
}

//Vec3D Vec3D::cross(const Vec3D &right) {
//	Vec3D temp;
//	temp.comp[0]=comp[1]*right.comp[2]-comp[2]*right.comp[1];
//	temp.comp[1]=-comp[0]*right.comp[2]+comp[2]*right.comp[0];
//	temp.comp[2]=comp[0]*right.comp[1]-comp[1]*right.comp[0];
//	return temp;
//}

VecND VecND::norm(void) {
	return (*this)/=fabs(*this);
}

VecND &VecND::operator= (const VecND &right) {
        comp=right.comp;
///	comp[0]=right.comp[0];
//	comp[1]=right.comp[1];
//	comp[2]=right.comp[2];
	return *this;
}

VecND &VecND::operator= (const double &right) {
        for(vdsz_type cnt=0; cnt!=comp.size();++cnt)
                 comp[cnt]=right;
	//comp[0]=right;
	//comp[1]=right;
	//comp[2]=right;
	return *this;
}

VecND &VecND::operator*= (const double &right) {
        for(vdsz_type cnt=0; cnt!=comp.size();++cnt)
                 comp[cnt]*=right;
	//comp[0]*=right;
//	comp[1]*=right;
//	comp[2]*=right;
	return *this;
}

VecND VecND::operator*(const double &right) {
	VecND temp(comp.size());
       for(vdsz_type cnt=0; cnt!=comp.size();++cnt) {
           temp.comp[cnt]=comp[cnt]*right;
        }
        	
//	Vec3D temp;
//	temp.comp[0]=comp[0]*right;
//	temp.comp[1]=comp[1]*right;
//	temp.comp[2]=comp[2]*right;
	return temp;
}

VecND &VecND::operator/= (const double &right) {
     for(vdsz_type cnt=0; cnt!=comp.size();++cnt) 
        comp[cnt]/=right;
	//comp[0]/=right;
	//comp[1]/=right;
	//comp[2]/=right;
	return *this;
}

VecND VecND::operator/ (const double &right) {
	VecND temp(comp.size());
       for(vdsz_type cnt=0; cnt!=comp.size();++cnt) 
             temp.comp[cnt]=comp[cnt]/right;
	//Vec3D temp;
//	temp.comp[0]=comp[0]/right;
//	temp.comp[1]=comp[1]/right;
//	temp.comp[2]=comp[2]/right;
	return temp;
}

VecND &VecND::operator+= (const double &right) {
       for(vdsz_type cnt=0; cnt!=comp.size();++cnt) 
          comp[cnt]+=right;
///	comp[0]+=right;
//	comp[1]+=right;
//	comp[2]+=right;
	return *this;
}

VecND &VecND::operator+= (const VecND &right) {
       for(vdsz_type cnt=0; cnt!=comp.size();++cnt) 
          comp[cnt]+=right.comp[cnt];
//	comp[0]+=right.comp[0];
//	comp[1]+=right.comp[1];
//	comp[2]+=right.comp[2];
	return *this;
}

VecND VecND::operator+ (const double &right) {
	VecND temp(comp.size());
       for(vdsz_type cnt=0; cnt!=comp.size();++cnt) 
           temp.comp[cnt]=comp[cnt]+right;
	//Vec3D temp;
	//temp.comp[0]=comp[0]+right;
	//temp.comp[1]=comp[1]+right;
	//temp.comp[2]=comp[2]+right;
	return temp;
}

VecND &VecND::operator-= (const double &right) {

       for(vdsz_type cnt=0; cnt!=comp.size();++cnt) 
         comp[cnt]-=right;
//	comp[0]-=right;
///	comp[1]-=right;
//	comp[2]-=right;
	return *this;
}

VecND &VecND::operator-= (const VecND &right) {

       for(vdsz_type cnt=0; cnt!=comp.size();++cnt) 
            comp[cnt]-=right.comp[cnt];
//	comp[0]-=right.comp[0];
//	comp[1]-=right.comp[1];
//	comp[2]-=right.comp[2];
	return *this;
}

VecND VecND::operator- (const double &right) {
	VecND temp(comp.size());
       for(vdsz_type cnt=0; cnt!=comp.size();++cnt) 
          temp.comp[cnt]=comp[cnt]-right;
	//VecND temp;
//	temp.comp[0]=comp[0]-right;
	//temp.comp[1]=comp[1]-right;
//	temp.comp[2]=comp[2]-right;
	return temp;
}

VecND operator*(const double &left, const VecND &right) {
	VecND temp(right.comp.size());
       for(vdsz_type cnt=0; cnt!=right.comp.size();++cnt) 
          temp.comp[cnt]=left*right.comp[cnt];
	//Vec3D temp;
	//temp.comp[0]=left*right.comp[0];
	//temp.comp[1]=left*right.comp[1];
//	temp.comp[2]=left*right.comp[2];
	return temp;
}

VecND operator/ (const double &left, const VecND &right) {
	VecND temp(right.comp.size());
       for(vdsz_type cnt=0; cnt!=right.comp.size();++cnt) 
          temp.comp[cnt]=left/right.comp[cnt];
//	Vec3D temp;
//	temp.comp[0]=left/right.comp[0];
//	temp.comp[1]=left/right.comp[1];
//	temp.comp[2]=left/right.comp[2];
	return temp;
}

VecND operator+ (const double &left, const VecND &right) {
	VecND temp(right.comp.size());
       for(vdsz_type cnt=0; cnt!=right.comp.size();++cnt) 
          temp.comp[cnt]=left+right.comp[cnt];
//	Vec3D temp;
//	temp.comp[0]=left+right.comp[0];
//	temp.comp[1]=left+right.comp[1];
//	temp.comp[2]=left+right.comp[2];
	return temp;
}

VecND operator+ (const VecND &left, const VecND &right) {
	VecND temp(right.comp.size());
       for(vdsz_type cnt=0; cnt!=right.comp.size();++cnt) 
          temp.comp[cnt]=left.comp[cnt]+right.comp[cnt];
//	Vec3D temp;
//	temp.comp[0]=left.comp[0]+right.comp[0];
//	temp.comp[1]=left.comp[1]+right.comp[1];
//	temp.comp[2]=left.comp[2]+right.comp[2];
	return temp;
}

VecND operator- (const double &left, const VecND &right) {
	VecND temp(right.comp.size());
       for(vdsz_type cnt=0; cnt!=right.comp.size();++cnt) 
            temp.comp[cnt]=left-right.comp[cnt];
//	Vec3D temp;
//	temp.comp[0]=left-right.comp[0];
//	temp.comp[1]=left-right.comp[1];
//	temp.comp[2]=left-right.comp[2];
	return temp;
}

VecND operator- (const VecND &left, const VecND &right) {
	VecND temp(right.comp.size());
       for(vdsz_type cnt=0; cnt!=right.comp.size();++cnt) 
           temp.comp[cnt]=left.comp[cnt]-right.comp[cnt];
//	Vec3D temp;
//	temp.comp[0]=left.comp[0]-right.comp[0];
//	temp.comp[1]=left.comp[1]-right.comp[1];
//	temp.comp[2]=left.comp[2]-right.comp[2];
	return temp;
}

bool VecND::operator== (const VecND &right) {
        return (comp==right.comp);
	//return (comp[0]==right.comp[0] && comp[1]==right.comp[1] && comp[2]==right.comp[2]);
}

bool VecND::operator!= (const VecND &right) {
        return (comp!=right.comp);
//	return (comp[0]!=right.comp[0] || comp[1]!=right.comp[1] || comp[2]!=right.comp[2]);
}

double &VecND::operator[] (int i) {return comp[i];}

ostream &operator<< (ostream &output,const VecND &right) {
	//output << "{" << right.comp[0] << "," << right.comp[1] << "," << right.comp[2] << "}";
	output << "{";
       for(vdsz_type cnt=0; cnt!=right.comp.size()-1;++cnt) 
         	output << right.comp[cnt]<<",";   
       output << right.comp[right.comp.size()-1]<<"}";  
	return output;
}
