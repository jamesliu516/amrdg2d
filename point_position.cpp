/* 射线法判断点q与多边形polygon的位置关系，要求polygon为简单多边形，顶点逆时针排列

                                                                    如果点在多边形内： 返回0

                                                                    如果点在多边形边上： 返回1

                                                                    如果点在多边形外： 返回2

*/

const double INFINITY = 1e10;

const double ESP = 1e-5;

const int MAX_N = 1000;



struct Point {

    double x, y;

};

struct LineSegment {

    Point pt1, pt2;

};

typedef vector<Point> Polygon;



// 计算叉乘 |P0P1| × |P0P2|

double Multiply(Point p1, Point p2, Point p0)

{

    return ( (p1.x - p0.x) * (p2.y - p0.y) - (p2.x - p0.x) * (p1.y - p0.y) );

}

// 判断线段是否包含点point

bool IsOnline(Point point, LineSegment line)

{

    return( ( fabs(Multiply(line.pt1, line.pt2, point)) < ESP ) &&

        ( ( point.x - line.pt1.x ) * ( point.x - line.pt2.x ) <= 0 ) &&

        ( ( point.y - line.pt1.y ) * ( point.y - line.pt2.y ) <= 0 ) );

}

// 判断线段相交

bool Intersect(LineSegment L1, LineSegment L2)

{

    return( (max(L1.pt1.x, L1.pt2.x) >= min(L2.pt1.x, L2.pt2.x)) &&

        (max(L2.pt1.x, L2.pt2.x) >= min(L1.pt1.x, L1.pt2.x)) &&

        (max(L1.pt1.y, L1.pt2.y) >= min(L2.pt1.y, L2.pt2.y)) &&

        (max(L2.pt1.y, L2.pt2.y) >= min(L1.pt1.y, L1.pt2.y)) &&

        (Multiply(L2.pt1, L1.pt2, L1.pt1) * Multiply(L1.pt2, L2.pt2, L1.pt1) >= 0) &&

        (Multiply(L1.pt1, L2.pt2, L2.pt1) * Multiply(L2.pt2, L1.pt2, L2.pt1) >= 0)

        );

}

// 判断点在多边形内

bool InPolygon(const Polygon& polygon, Point point)

{

    int n = polygon.size();

    int count = 0;

    LineSegment line;

    line.pt1 = point;

    line.pt2.y = point.y;

    line.pt2.x = - INFINITY;



    for( int i = 0; i < n; i++ ) {

        // 得到多边形的一条边

        LineSegment side;

        side.pt1 = polygon[i];

        side.pt2 = polygon[(i + 1) % n];



        if( IsOnline(point, side) ) {

            return 1 ;

        }



        // 如果side平行x轴则不作考虑

        if( fabs(side.pt1.y - side.pt2.y) < ESP ) {

            continue;

        }



        if( IsOnline(side.pt1, line) ) {

            if( side.pt1.y > side.pt2.y ) count++;

        } else if( IsOnline(side.pt2, line) ) {

            if( side.pt2.y > side.pt1.y ) count++;

        } else if( Intersect(line, side) ) {

            count++;

        }

    }



    if ( count % 2 == 1 ) {return 0;}

    else { return 2;}

}

}
