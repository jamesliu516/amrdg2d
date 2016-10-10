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
#ifndef VEC4D_H
#define VEC4D_H

#include <iostream>
using namespace std;

class Vec4D {
public:
	double comp[4];  //点的坐标，comp[0]:x,comp[1]:y, comp[2]:z
	Vec4D(double x=0., double y=0., double z=0., double r=0.0);
	double dot(const Vec4D &right);
//	Vec4D cross(const Vec4D &right);
	Vec4D norm(void);//单位向量
	Vec4D &operator= (const Vec4D &);
	Vec4D &operator= (const double &);
	Vec4D &operator*= (const double &);
	Vec4D operator*(const double &);
	Vec4D &operator/= (const double &);
	Vec4D operator/ (const double &);
	Vec4D &operator+= (const double &);
	Vec4D &operator+= (const Vec4D &);
	Vec4D operator+ (const double &);
	Vec4D &operator-= (const double &);
	Vec4D &operator-= (const Vec4D &);
	Vec4D operator- (const double &);
	bool operator== (const Vec4D &);
	bool operator!= (const Vec4D &);
	double &operator[] (int i); //使用下标方式提取comp[0-3]
};

double fabs(const Vec4D vec);
Vec4D operator*(const double &left, const Vec4D &right);
Vec4D operator/ (const double &left, const Vec4D &right);
Vec4D operator+ (const double &left, const Vec4D &right);
Vec4D operator- (const double &left, const Vec4D &right);

Vec4D operator+ (const Vec4D &left, const Vec4D &right);
Vec4D operator- (const Vec4D &left, const Vec4D &right);

ostream &operator<< (ostream &output,const Vec4D &right);

#endif
