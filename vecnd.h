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
#ifndef VECND_H
#define VECND_H

#include <iostream>
#include<vector>

using namespace std;

typedef vector<double>::size_type vdsz_type;

class VecND {
public:
	//double comp[3];  //点的坐标，comp[0]:x,comp[1]:y, comp[2]:z
	vector<double> comp;
//	Vec3D(double x=0., double y=0., double z=0.);
       VecND(int n=1);
       
	double dot(const VecND &right);
	//Vec3D cross(const Vec3D &right);
	VecND norm(void);
	VecND &operator= (const VecND &);
	VecND &operator= (const double &);
	VecND &operator*= (const double &);
	VecND operator*(const double &);
	VecND &operator/= (const double &);
	VecND operator/ (const double &);
	VecND &operator+= (const double &);
	VecND &operator+= (const VecND &);
	VecND operator+ (const double &);
	VecND &operator-= (const double &);
	VecND &operator-= (const VecND &);
	VecND operator- (const double &);
	bool operator== (const VecND &);
	bool operator!= (const VecND &);
	double &operator[] (int i); //使用下标方式提取comp[0-2]
	void resize(int n, double r=0.0);
	int size();
};

double fabs(const VecND &vec);
VecND operator*(const double &left, const VecND &right);
VecND operator/ (const double &left, const VecND &right);
VecND operator+ (const double &left, const VecND &right);
VecND operator- (const double &left, const VecND &right);

VecND operator+ (const VecND &left, const VecND &right);
VecND operator- (const VecND &left, const VecND &right);

ostream &operator<< (ostream &output,const VecND &right);

#endif
