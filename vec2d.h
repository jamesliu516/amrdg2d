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
#ifndef VEC2D_H
#define VEC2D_H

#include <iostream>
using namespace std;

class Vec2D {
public:
	double comp[2];  //点的坐标，comp[0]:x,comp[1]:y, comp[2]:z
	//Vec3D(double x=0., double y=0., double z=0.);
	Vec2D(double x=0., double y=0.);
	double dot(const Vec2D &right);
//	Vec2D cross(const Vec2D &right);
	Vec2D norm(void);
	Vec2D &operator= (const Vec2D &);
	Vec2D &operator= (const double &);
	Vec2D &operator*= (const double &);
	Vec2D operator*(const double &);
	Vec2D &operator/= (const double &);
	Vec2D operator/ (const double &);
	Vec2D &operator+= (const double &);
	Vec2D &operator+= (const Vec2D &);
	Vec2D operator+ (const double &);
	Vec2D &operator-= (const double &);
	Vec2D &operator-= (const Vec2D &);
	Vec2D operator- (const double &);
	bool operator== (const Vec2D &);
	bool operator!= (const Vec2D &);	
	double &operator[] (int i); //使用下标方式提取comp[0-2]
};

double fabs(const Vec2D vec);
Vec2D operator*(const double &left, const Vec2D &right);
Vec2D operator/ (const double &left, const Vec2D &right);
Vec2D operator+ (const double &left, const Vec2D &right);
Vec2D operator- (const double &left, const Vec2D &right);

Vec2D operator+ (const Vec2D &left, const Vec2D &right);
Vec2D operator- (const Vec2D &left, const Vec2D &right);

ostream &operator<< (ostream &output,const Vec2D &right);

#endif
