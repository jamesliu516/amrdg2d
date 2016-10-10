/************************************************************************
	
	Copyright 2007-2009 Emre Sozer & Patrick Clark Trizila

	Contact: emresozer@freecfd.com , ptrizila@freecfd.com

	This file is a part of Free CFD

	Free CFD is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
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
#include "vec2d.h"

Vec2D::Vec2D(double x, double y) {
	comp[0]=x;
	comp[1]=y;
	//comp[2]=z;
}

double Vec2D::dot(const Vec2D &right) {
	//return (comp[0]*right.comp[0]+comp[1]*right.comp[1]+comp[2]*right.comp[2]);
	return (comp[0]*right.comp[0]+comp[1]*right.comp[1]);
}

double fabs(const Vec2D vec) {
	//return sqrt(vec.comp[0]*vec.comp[0]+vec.comp[1]*vec.comp[1]+vec.comp[2]*vec.comp[2]);
	return sqrt(vec.comp[0]*vec.comp[0]+vec.comp[1]*vec.comp[1]);
}
/*
Vec3D Vec3D::cross(const Vec3D &right) {
	Vec3D temp;
	temp.comp[0]=comp[1]*right.comp[2]-comp[2]*right.comp[1];
	temp.comp[1]=-comp[0]*right.comp[2]+comp[2]*right.comp[0];
	temp.comp[2]=comp[0]*right.comp[1]-comp[1]*right.comp[0];
	return temp;
}
*/

Vec2D Vec2D::norm(void) {
	return (*this)/=fabs(*this);
}

Vec2D &Vec2D::operator= (const Vec2D &right) {
	comp[0]=right.comp[0];
	comp[1]=right.comp[1];
//	comp[2]=right.comp[2];
	return *this;
}

Vec2D &Vec2D::operator= (const double &right) {
	comp[0]=right;
	comp[1]=right;
//	comp[2]=right;
	return *this;
}

Vec2D &Vec2D::operator*= (const double &right) {
	comp[0]*=right;
	comp[1]*=right;
//	comp[2]*=right;
	return *this;
}

Vec2D Vec2D::operator*(const double &right) {
	Vec2D temp;
	temp.comp[0]=comp[0]*right;
	temp.comp[1]=comp[1]*right;
//	temp.comp[2]=comp[2]*right;
	return temp;
}

Vec2D &Vec2D::operator/= (const double &right) {
	comp[0]/=right;
	comp[1]/=right;
	//comp[2]/=right;
	return *this;
}

Vec2D Vec2D::operator/ (const double &right) {
	Vec2D temp;
	temp.comp[0]=comp[0]/right;
	temp.comp[1]=comp[1]/right;
//	temp.comp[2]=comp[2]/right;
	return temp;
}

Vec2D &Vec2D::operator+= (const double &right) {
	comp[0]+=right;
	comp[1]+=right;
//	comp[2]+=right;
	return *this;
}

Vec2D &Vec2D::operator+= (const Vec2D &right) {
	comp[0]+=right.comp[0];
	comp[1]+=right.comp[1];
//	comp[2]+=right.comp[2];
	return *this;
}

Vec2D Vec2D::operator+ (const double &right) {
	Vec2D temp;
	temp.comp[0]=comp[0]+right;
	temp.comp[1]=comp[1]+right;
//	temp.comp[2]=comp[2]+right;
	return temp;
}

Vec2D &Vec2D::operator-= (const double &right) {
	comp[0]-=right;
	comp[1]-=right;
//	comp[2]-=right;
	return *this;
}

Vec2D &Vec2D::operator-= (const Vec2D &right) {
	comp[0]-=right.comp[0];
	comp[1]-=right.comp[1];
	//comp[2]-=right.comp[2];
	return *this;
}

Vec2D Vec2D::operator- (const double &right) {
	Vec2D temp;
	temp.comp[0]=comp[0]-right;
	temp.comp[1]=comp[1]-right;
	//temp.comp[2]=comp[2]-right;
	return temp;
}

Vec2D operator*(const double &left, const Vec2D &right) {
	Vec2D temp;
	temp.comp[0]=left*right.comp[0];
	temp.comp[1]=left*right.comp[1];
//	temp.comp[2]=left*right.comp[2];
	return temp;
}

Vec2D operator/ (const double &left, const Vec2D &right) {
	Vec2D temp;
	temp.comp[0]=left/right.comp[0];
	temp.comp[1]=left/right.comp[1];
//	temp.comp[2]=left/right.comp[2];
	return temp;
}

Vec2D operator+ (const double &left, const Vec2D &right) {
	Vec2D temp;
	temp.comp[0]=left+right.comp[0];
	temp.comp[1]=left+right.comp[1];
//	temp.comp[2]=left+right.comp[2];
	return temp;
}

Vec2D operator+ (const Vec2D &left, const Vec2D &right) {
	Vec2D temp;
	temp.comp[0]=left.comp[0]+right.comp[0];
	temp.comp[1]=left.comp[1]+right.comp[1];
//	temp.comp[2]=left.comp[2]+right.comp[2];
	return temp;
}

Vec2D operator- (const double &left, const Vec2D &right) {
	Vec2D temp;
	temp.comp[0]=left-right.comp[0];
	temp.comp[1]=left-right.comp[1];
//	temp.comp[2]=left-right.comp[2];
	return temp;
}

Vec2D operator- (const Vec2D &left, const Vec2D &right) {
	Vec2D temp;
	temp.comp[0]=left.comp[0]-right.comp[0];
	temp.comp[1]=left.comp[1]-right.comp[1];
//	temp.comp[2]=left.comp[2]-right.comp[2];
	return temp;
}

bool Vec2D::operator== (const Vec2D &right) {
	return (comp[0]==right.comp[0] && comp[1]==right.comp[1]);
}

bool Vec2D::operator!= (const Vec2D &right) {
	return (comp[0]!=right.comp[0] || comp[1]!=right.comp[1]);
}

double &Vec2D::operator[] (int i) {return comp[i];}

ostream &operator<< (ostream &output,const Vec2D &right) {
	output << "{" << right.comp[0] << "," << right.comp[1]  << "}";
	return output;
}
