/************************************************************************
	
	Copyright 2007-2009 Emre Sozer & Patrick Clark Trizila

	Contact: emresozer@freecfd.com , ptrizila@freecfd.com

	This file is a part of Free CFD

	Free CFD is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 4 of the License, or
    any later version.

    Free CFD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    For a copy of the GNU General Public License,
    see <http://www.gnu.org/licenses/>.

*************************************************************************/
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;
#include "vec4d.h"

Vec4D::Vec4D(double x, double y, double z, double r) {
	comp[0]=x;
	comp[1]=y;
	comp[2]=z;
       comp[3]=r;
}

double Vec4D::dot(const Vec4D &right) {
	return (comp[0]*right.comp[0]+comp[1]*right.comp[1]+comp[2]*right.comp[2]+comp[3]*right.comp[3]);
}

double fabs(const Vec4D vec) {
	return sqrt(vec.comp[0]*vec.comp[0]+vec.comp[1]*vec.comp[1]+vec.comp[2]*vec.comp[2]+vec.comp[3]*vec.comp[3]);
}
/*
Vec4D Vec4D::cross(const Vec4D &right) {
	Vec4D temp;
	temp.comp[0]=comp[1]*right.comp[2]-comp[2]*right.comp[1];
	temp.comp[1]=-comp[0]*right.comp[2]+comp[2]*right.comp[0];
	temp.comp[2]=comp[0]*right.comp[1]-comp[1]*right.comp[0];
	return temp;
}*/

Vec4D Vec4D::norm(void) {
	return (*this)/=fabs(*this);
}

Vec4D &Vec4D::operator= (const Vec4D &right) {
	comp[0]=right.comp[0];
	comp[1]=right.comp[1];
	comp[2]=right.comp[2];
	comp[3]=right.comp[3];	
	return *this;
}

Vec4D &Vec4D::operator= (const double &right) {
	comp[0]=right;
	comp[1]=right;
	comp[2]=right;
	comp[3]=right;
	return *this;
}

Vec4D &Vec4D::operator*= (const double &right) {
	comp[0]*=right;
	comp[1]*=right;
	comp[2]*=right;
	comp[3]*=right;	
	
	return *this;
}

Vec4D Vec4D::operator*(const double &right) {
	Vec4D temp;
	temp.comp[0]=comp[0]*right;
	temp.comp[1]=comp[1]*right;
	temp.comp[2]=comp[2]*right;
	temp.comp[3]=comp[3]*right;	
	return temp;
}

Vec4D &Vec4D::operator/= (const double &right) {
	comp[0]/=right;
	comp[1]/=right;
	comp[2]/=right;
	comp[3]/=right;	
	return *this;
}

Vec4D Vec4D::operator/ (const double &right) {
	Vec4D temp;
	temp.comp[0]=comp[0]/right;
	temp.comp[1]=comp[1]/right;
	temp.comp[2]=comp[2]/right;
	temp.comp[3]=comp[3]/right;	
	
	return temp;
}

Vec4D &Vec4D::operator+= (const double &right) {
	comp[0]+=right;
	comp[1]+=right;
	comp[2]+=right;
	comp[3]+=right;	
	return *this;
}

Vec4D &Vec4D::operator+= (const Vec4D &right) {
	comp[0]+=right.comp[0];
	comp[1]+=right.comp[1];
	comp[2]+=right.comp[2];
	comp[3]+=right.comp[3];	
	return *this;
}

Vec4D Vec4D::operator+ (const double &right) {
	Vec4D temp;
	temp.comp[0]=comp[0]+right;
	temp.comp[1]=comp[1]+right;
	temp.comp[2]=comp[2]+right;
	temp.comp[3]=comp[3]+right;	
	return temp;
}

Vec4D &Vec4D::operator-= (const double &right) {
	comp[0]-=right;
	comp[1]-=right;
	comp[2]-=right;
	comp[3]-=right;	
	return *this;
}

Vec4D &Vec4D::operator-= (const Vec4D &right) {
	comp[0]-=right.comp[0];
	comp[1]-=right.comp[1];
	comp[2]-=right.comp[2];
	comp[3]-=right.comp[3];	
	return *this;
}

Vec4D Vec4D::operator- (const double &right) {
	Vec4D temp;
	temp.comp[0]=comp[0]-right;
	temp.comp[1]=comp[1]-right;
	temp.comp[2]=comp[2]-right;
       temp.comp[3]=comp[3]-right;
	return temp;
}

Vec4D operator*(const double &left, const Vec4D &right) {
	Vec4D temp;
	temp.comp[0]=left*right.comp[0];
	temp.comp[1]=left*right.comp[1];
	temp.comp[2]=left*right.comp[2];
	temp.comp[3]=left*right.comp[3];
	return temp;
}

Vec4D operator/ (const double &left, const Vec4D &right) {
	Vec4D temp;
	temp.comp[0]=left/right.comp[0];
	temp.comp[1]=left/right.comp[1];
	temp.comp[2]=left/right.comp[2];
	temp.comp[3]=left/right.comp[3];
	return temp;
}

Vec4D operator+ (const double &left, const Vec4D &right) {
	Vec4D temp;
	temp.comp[0]=left+right.comp[0];
	temp.comp[1]=left+right.comp[1];
	temp.comp[2]=left+right.comp[2];
	temp.comp[3]=left+right.comp[3];
	return temp;
}

Vec4D operator+ (const Vec4D &left, const Vec4D &right) {
	Vec4D temp;
	temp.comp[0]=left.comp[0]+right.comp[0];
	temp.comp[1]=left.comp[1]+right.comp[1];
	temp.comp[2]=left.comp[2]+right.comp[2];
	temp.comp[3]=left.comp[3]+right.comp[3];
	return temp;
}

Vec4D operator- (const double &left, const Vec4D &right) {
	Vec4D temp;
	temp.comp[0]=left-right.comp[0];
	temp.comp[1]=left-right.comp[1];
	temp.comp[2]=left-right.comp[2];
	temp.comp[3]=left-right.comp[3];
	return temp;
}

Vec4D operator- (const Vec4D &left, const Vec4D &right) {
	Vec4D temp;
	temp.comp[0]=left.comp[0]-right.comp[0];
	temp.comp[1]=left.comp[1]-right.comp[1];
	temp.comp[2]=left.comp[2]-right.comp[2];
	temp.comp[3]=left.comp[3]-right.comp[3];
	return temp;
}

bool Vec4D::operator== (const Vec4D &right) {
	return (comp[0]==right.comp[0] && comp[1]==right.comp[1] && comp[2]==right.comp[2]&& comp[3]==right.comp[3]);
}

bool Vec4D::operator!= (const Vec4D &right) {
	return (comp[0]!=right.comp[0] || comp[1]!=right.comp[1] || comp[2]!=right.comp[2]|| comp[3]!=right.comp[3]);
}

double &Vec4D::operator[] (int i) {return comp[i];}

ostream &operator<< (ostream &output,const Vec4D &right) {
	output << "{" << right.comp[0] << "," << right.comp[1] << "," << right.comp[2] << "," << right.comp[3]<< "}";
	return output;
}
