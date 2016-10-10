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
#ifndef VEC3D_H
#define VEC3D_H

#include <iostream>
using namespace std;

class Vec3D {
public:
	double comp[3];  //点的坐标，comp[0]:x,comp[1]:y, comp[2]:z
	Vec3D(double x=0., double y=0., double z=0.);
	double dot(const Vec3D &right);
	Vec3D cross(const Vec3D &right);
	Vec3D norm(void);
	Vec3D &operator= (const Vec3D &);
	Vec3D &operator= (const double &);
	Vec3D &operator*= (const double &);
	Vec3D operator*(const double &);
	Vec3D &operator/= (const double &);
	Vec3D operator/ (const double &);
	Vec3D &operator+= (const double &);
	Vec3D &operator+= (const Vec3D &);
	Vec3D operator+ (const double &);
	Vec3D &operator-= (const double &);
	Vec3D &operator-= (const Vec3D &);
	Vec3D operator- (const double &);
	bool operator== (const Vec3D &);
	bool operator!= (const Vec3D &);
	double &operator[] (int i); //使用下标方式提取comp[0-2]
};

double fabs(const Vec3D vec);
Vec3D operator*(const double &left, const Vec3D &right);
Vec3D operator/ (const double &left, const Vec3D &right);
Vec3D operator+ (const double &left, const Vec3D &right);
Vec3D operator- (const double &left, const Vec3D &right);

Vec3D operator+ (const Vec3D &left, const Vec3D &right);
Vec3D operator- (const Vec3D &left, const Vec3D &right);

ostream &operator<< (ostream &output,const Vec3D &right);

#endif
