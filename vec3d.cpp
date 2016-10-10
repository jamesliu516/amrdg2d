/************************************************************************
	
	Copyright 2007-2009 Emre Sozer & Patrick Clark Trizila

	Contact: emresozer@freecfd.com , ptrizila@freecfd.com

	This file is a part of Free CFD

	Free CFD is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
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
#include "vec3d.h"

Vec3D::Vec3D(double x, double y, double z) {
	comp[0]=x;
	comp[1]=y;
	comp[2]=z;
}

double Vec3D::dot(const Vec3D &right) {
	return (comp[0]*right.comp[0]+comp[1]*right.comp[1]+comp[2]*right.comp[2]);
}

double fabs(const Vec3D vec) {
	return sqrt(vec.comp[0]*vec.comp[0]+vec.comp[1]*vec.comp[1]+vec.comp[2]*vec.comp[2]);
}

Vec3D Vec3D::cross(const Vec3D &right) {
	Vec3D temp;
	temp.comp[0]=comp[1]*right.comp[2]-comp[2]*right.comp[1];
	temp.comp[1]=-comp[0]*right.comp[2]+comp[2]*right.comp[0];
	temp.comp[2]=comp[0]*right.comp[1]-comp[1]*right.comp[0];
	return temp;
}

Vec3D Vec3D::norm(void) {
	return (*this)/=fabs(*this);
}

Vec3D &Vec3D::operator= (const Vec3D &right) {
	comp[0]=right.comp[0];
	comp[1]=right.comp[1];
	comp[2]=right.comp[2];
	return *this;
}

Vec3D &Vec3D::operator= (const double &right) {
	comp[0]=right;
	comp[1]=right;
	comp[2]=right;
	return *this;
}

Vec3D &Vec3D::operator*= (const double &right) {
	comp[0]*=right;
	comp[1]*=right;
	comp[2]*=right;
	return *this;
}

Vec3D Vec3D::operator*(const double &right) {
	Vec3D temp;
	temp.comp[0]=comp[0]*right;
	temp.comp[1]=comp[1]*right;
	temp.comp[2]=comp[2]*right;
	return temp;
}

Vec3D &Vec3D::operator/= (const double &right) {
	comp[0]/=right;
	comp[1]/=right;
	comp[2]/=right;
	return *this;
}

Vec3D Vec3D::operator/ (const double &right) {
	Vec3D temp;
	temp.comp[0]=comp[0]/right;
	temp.comp[1]=comp[1]/right;
	temp.comp[2]=comp[2]/right;
	return temp;
}

Vec3D &Vec3D::operator+= (const double &right) {
	comp[0]+=right;
	comp[1]+=right;
	comp[2]+=right;
	return *this;
}

Vec3D &Vec3D::operator+= (const Vec3D &right) {
	comp[0]+=right.comp[0];
	comp[1]+=right.comp[1];
	comp[2]+=right.comp[2];
	return *this;
}

Vec3D Vec3D::operator+ (const double &right) {
	Vec3D temp;
	temp.comp[0]=comp[0]+right;
	temp.comp[1]=comp[1]+right;
	temp.comp[2]=comp[2]+right;
	return temp;
}

Vec3D &Vec3D::operator-= (const double &right) {
	comp[0]-=right;
	comp[1]-=right;
	comp[2]-=right;
	return *this;
}

Vec3D &Vec3D::operator-= (const Vec3D &right) {
	comp[0]-=right.comp[0];
	comp[1]-=right.comp[1];
	comp[2]-=right.comp[2];
	return *this;
}

Vec3D Vec3D::operator- (const double &right) {
	Vec3D temp;
	temp.comp[0]=comp[0]-right;
	temp.comp[1]=comp[1]-right;
	temp.comp[2]=comp[2]-right;
	return temp;
}

Vec3D operator*(const double &left, const Vec3D &right) {
	Vec3D temp;
	temp.comp[0]=left*right.comp[0];
	temp.comp[1]=left*right.comp[1];
	temp.comp[2]=left*right.comp[2];
	return temp;
}

Vec3D operator/ (const double &left, const Vec3D &right) {
	Vec3D temp;
	temp.comp[0]=left/right.comp[0];
	temp.comp[1]=left/right.comp[1];
	temp.comp[2]=left/right.comp[2];
	return temp;
}

Vec3D operator+ (const double &left, const Vec3D &right) {
	Vec3D temp;
	temp.comp[0]=left+right.comp[0];
	temp.comp[1]=left+right.comp[1];
	temp.comp[2]=left+right.comp[2];
	return temp;
}

Vec3D operator+ (const Vec3D &left, const Vec3D &right) {
	Vec3D temp;
	temp.comp[0]=left.comp[0]+right.comp[0];
	temp.comp[1]=left.comp[1]+right.comp[1];
	temp.comp[2]=left.comp[2]+right.comp[2];
	return temp;
}

Vec3D operator- (const double &left, const Vec3D &right) {
	Vec3D temp;
	temp.comp[0]=left-right.comp[0];
	temp.comp[1]=left-right.comp[1];
	temp.comp[2]=left-right.comp[2];
	return temp;
}

Vec3D operator- (const Vec3D &left, const Vec3D &right) {
	Vec3D temp;
	temp.comp[0]=left.comp[0]-right.comp[0];
	temp.comp[1]=left.comp[1]-right.comp[1];
	temp.comp[2]=left.comp[2]-right.comp[2];
	return temp;
}

bool Vec3D::operator== (const Vec3D &right) {
	return (comp[0]==right.comp[0] && comp[1]==right.comp[1] && comp[2]==right.comp[2]);
}

bool Vec3D::operator!= (const Vec3D &right) {
	return (comp[0]!=right.comp[0] || comp[1]!=right.comp[1] || comp[2]!=right.comp[2]);
}

double &Vec3D::operator[] (int i) {return comp[i];}

ostream &operator<< (ostream &output,const Vec3D &right) {
	output << "{" << right.comp[0] << "," << right.comp[1] << "," << right.comp[2] << "}";
	return output;
}
