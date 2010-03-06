#include"ellipse-mat.h"


//****************************************************************************
//***  CLASSES TO CONSTRUCT ELLIPSES AND ADDITIONAL RELATED FUNCTIONS  *******
//*** Similar to the classe CEllipse but with (x,y,z) changed in (0,x,y)   ***
//*** and alpha=0      *******************************************************
//***  'ratio'=b/a for the ellipse equation
//***     x^2/a^2+y^2/b^2=1
//***  with a>b

//****************************************************************************

CEllipse::CEllipse(double x,double y,double r,double ratio,double _beta, int c, int flags) {setup();};

void CEllipse::setup(){

b=e*radius;     //length of minor axis   

rotat_mat(2,2)= 1;//z component

scale_mat(2,2)= 1;

ellip_mat(2,2)= 1;

//Elements of the rotational matrix

  rotat_mat(0,0)= cos(beta);
  rotat_mat(0,1)= -sin(beta);
  rotat_mat(1,0)= sin(beta);
  rotat_mat(1,1)= cos(beta);

  //Elements of the scaling matrix

  scale_mat(0,0)=1.0/(radius*radius);
  scale_mat(1,1)=1.0/(b*b);

  ellip_mat=rotat_mat*scale_mat*~rotat_mat;

Matrix test(3,3);
cerr<< "mat "<<test <<endl;
}


void CEllipse::set_ellipse(double x, double y , double r, double _e, double _beta, int c, int flags){
		set(x,y,r,c,flags);	
		e=_e;
		beta=_beta;
		setup();
		}	
/*
***********************************************************************
****** FUNCTION TO PUT ONE ELLIPSE 'E' INTO CONTACT TO ANOTHER  ******
****** ONE 'E0' THROUGH A TRANSLATION         *************************
*** INPUT: 'E' is the *new* ellipse, the one to translate    *********
***        'E0' is the old ellipse that remains fixed         *********
*** OUTPUT: Vector 'vector_for_rescaling':                       *****
***         First component '.x' is the Lagrange multiplier       *****
***         Second component '.y' is the function that rescales   *****
***            the new ellipse 'E'.                    ***************
***         Third component '.z' is '1' if the nonlinear Newton-Raph.**
***            procedure converged, and '0' otherwise.  ***************
***********************************************************************
*/

/*
CVector rescale_ellipse_to_touch(const CEllipse &E, const CEllipse &E0){
//here we should if E will touch E0 or not;

Matrix R(3,1), R0(3,1), R12(3,1);
R12(0,0)=E0.x(0)-E.x(0);
R12(1,0)=E0.x(1)-E.x(1);
R12(2,0)=0;
R(0,0)=E.x(0);
R(1,0)=E.x(1);
R(2,0)=0;
R0(0,0)=E0.x(0);
R0(1,0)=E0.x(1);
R0(2,0)=0;


// Begin of iteration for lambda
Matrix D(3,3), invD(3,3), F(3,3), A(3,3), B(3,3);
Matrix invE(3,3), invE0(3,3);
invE=!E.ellip_mat;
invE0=!E0.ellip_mat;
int itermax=10;
double lambda=0.5;
for(int loop=1;loop<=itermax;loop++){

	D=invE0+lambda*invE;	
	invD=!D;
//these two lines can be optimized
	double d=lambda*lambda*invE.Det()+lambda*((invE+invE0).Det()-invE.Det()-invE0.Det())+invE0.Det();
	double dd=2*lambda*invE.Det()+((invE+invE0).Det()-invE.Det()-invE0.Det());
	
	F=invE*invE0*invD;
	A=invD*invE0*invD;
	
	B=2.0/d*(F-dd*A);



	double dx= ((~R)*A*R)(0,0)/((~R)*B*R)(0,0);

	double threshold=1e-5;
	Matrix E12(3,3);
	Matrix rc(3,1);
	Matrix rE(3,1);
	if(fabs(dx) < threshold){
		//CVector rc;
		E12=!(E.ellip_mat-(lambda*E0.ellip_mat));
		rc=R0-E12*(E.ellip_mat*R12);
		double s=lambda*sqrt((~R*(E0.ellip_mat*E12*E.ellip_mat*E12*E0.ellip_mat)*R)(0,0));
		rE=(1-1/s)*rc+1/s*R;
		return CVector(rE(0,0),rE(1,0));
		}

	lambda -= dx;
	}   

  return E.center;   //NON Convergence.

}
*/
