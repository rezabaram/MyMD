#ifndef ELLIPSE_H
#define ELLIPSE_H

#include"particle.h"
#include"matrix.h"

using namespace std;
using namespace math;

typedef matrix<double> Matrix;

class CEllipse {
	public:
	CEllipse(vec _x, vec _r, quat _q){set_ellipse(_x, _r, _q)};
	CEllipse():x(1.0), r(0), q(0){};
	void setup();
	void set_ellipse(vec _x, vec _r, quat _q);

	Matrix rotat_mat;
	Matrix scale_mat;
	Matrix ellip_mat;
	vec x;//center of mass
	quat q; //quaternion
	private:
	vec r;// radius and b are the major and minor axis of the ellipse
};

//CVector rescale_ellipse_to_touch(const CEllipse &E1, const CEllipse &E0);
//CVector touch_Point(const CEllipse &EE1, const CEllipse &EE0, const double &lambda);
//CVector NewCentroid_of_(const CEllipse &E0, const double RescalFact, const CVector TouchPoint);

//CVector rotate_ellipse_to_touch(const CEllipse &E1, const CEllipse &E0, const CVector &TouchPoint);
#endif
