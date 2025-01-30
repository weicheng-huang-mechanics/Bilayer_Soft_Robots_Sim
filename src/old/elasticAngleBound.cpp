#include "elasticAngleBound.h"
#include <iostream>

elasticAngleBound::elasticAngleBound(elasticPlate &m_plate, timeStepper &m_stepper)
{
	plate = &m_plate;
    stepper = &m_stepper;

	EA = plate->EA;
}

elasticAngleBound::~elasticAngleBound()
{
	;
}

void elasticAngleBound::computeFa()
{
	for(int k=0; k < plate->v_coupleAngle.size(); k++)
    {
    	int ind1 = plate->v_coupleAngle[k].nv_1;
        int ind2 = plate->v_coupleAngle[k].nv_2;
        int ind3 = plate->v_coupleAngle[k].nv_3;
        int ind4 = plate->v_coupleAngle[k].nv_4;

        Vector3d x1 = plate->getVertex(ind1);
        Vector3d x2 = plate->getVertex(ind2);
        Vector3d x3 = plate->getVertex(ind3);
        Vector3d x4 = plate->getVertex(ind4);

        VectorXd force = - EA * computeAngleForce(x1(0), x1(1), x1(2), x2(0), x2(1), x2(2), x3(0), x3(1), x3(2), x4(0), x4(1), x4(2));

        VectorXi localDOF = VectorXi::Zero(12);

        localDOF = plate->v_coupleAngle[k].arrayNum;

        for(int i = 0; i < 12; i++)
		{
			int ind = localDOF(i);
			stepper->addForce(ind, - force[i]);
		}
    }
}

void elasticAngleBound::computeJa()
{
	for(int k=0; k < plate->v_coupleAngle.size(); k++)
    {
    	int ind1 = plate->v_coupleAngle[k].nv_1;
        int ind2 = plate->v_coupleAngle[k].nv_2;
        int ind3 = plate->v_coupleAngle[k].nv_3;
        int ind4 = plate->v_coupleAngle[k].nv_4;

        Vector3d x1 = plate->getVertex(ind1);
        Vector3d x2 = plate->getVertex(ind2);
        Vector3d x3 = plate->getVertex(ind3);
        Vector3d x4 = plate->getVertex(ind4);

        MatrixXd jacob = EA * computeAngleJacobian(x1(0), x1(1), x1(2), x2(0), x2(1), x2(2), x3(0), x3(1), x3(2), x4(0), x4(1), x4(2));

        VectorXi localDOF = VectorXi::Zero(12);

        localDOF = plate->v_coupleAngle[k].arrayNum;

        for(int i = 0; i < 12; i++)
		{
			for(int j = 0; j < 12; j++)
			{
				int ind11 = localDOF(i);
				int ind22 = localDOF(j);
				
				stepper->addJacobian(ind11, ind22, jacob(i,j));
			}
		}
    }
}

void elasticAngleBound::setFirstJacobian()
{
	for(int k=0; k < plate->v_coupleAngle.size(); k++)
    {
        VectorXi localDOF = VectorXi::Zero(12);

        localDOF = plate->v_coupleAngle[k].arrayNum;

        for (int i = 0; i < 12; i++)
        {
            for (int j = 0; j < 12; j++)
            {
                stepper->addJacobian(localDOF(i), localDOF(j), 1);
            }
        }
    }
}

VectorXd elasticAngleBound::computeAngleForce(double x1, double y1, double z1, double x2, double y2, double z2, 
        double x3, double y3, double z3, double xe, double ye, double ze)
{
    VectorXd vecResult;

    vecResult = ListVec((1.*(x2 - xe)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + 
         z1*z2 - z2*z3 - z1*ze + z3*ze))/
     ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
       (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
    (1.*(x1 - x3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
         y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
     (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
       (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
   (1.*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + 
         z1*z2 - z2*z3 - z1*ze + z3*ze))/
     ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
       (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
    (1.*(y1 - y3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
         y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
     (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
       (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
   (1.*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + 
         z1*z2 - z2*z3 - z1*ze + z3*ze))/
     ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
       (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
    (1.*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
         y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
     (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
       (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
   (1.*(x1 - x3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + 
         z1*z2 - z2*z3 - z1*ze + z3*ze))/
     ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
       (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
    (1.*(x2 - xe)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
         y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
     ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
       pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
   (1.*(y1 - y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + 
         z1*z2 - z2*z3 - z1*ze + z3*ze))/
     ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
       (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
    (1.*(y2 - ye)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
         y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
     ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
       pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
   (1.*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + 
         z1*z2 - z2*z3 - z1*ze + z3*ze))/
     ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
       (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
    (1.*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
         y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
     ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
       pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
   (1.*(-x2 + xe)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + 
         z1*z2 - z2*z3 - z1*ze + z3*ze))/
     ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
       (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
    (1.*(x1 - x3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
         y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
     (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
       (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
   (1.*(-y2 + ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + 
         z1*z2 - z2*z3 - z1*ze + z3*ze))/
     ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
       (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
    (1.*(y1 - y3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
         y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
     (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
       (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
   (1.*(-z2 + ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + 
         z1*z2 - z2*z3 - z1*ze + z3*ze))/
     ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
       (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
    (1.*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
         y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
     (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
       (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
   (1.*(-x1 + x3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + 
         z1*z2 - z2*z3 - z1*ze + z3*ze))/
     ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
       (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
    (1.*(x2 - xe)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
         y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
     ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
       pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
   (1.*(-y1 + y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + 
         z1*z2 - z2*z3 - z1*ze + z3*ze))/
     ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
       (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
    (1.*(y2 - ye)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
         y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
     ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
       pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
   (1.*(-z1 + z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + 
         z1*z2 - z2*z3 - z1*ze + z3*ze))/
     ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
       (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
    (1.*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
         y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
     ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
       pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)));

    return vecResult;
}

MatrixXd elasticAngleBound::computeAngleJacobian(double x1, double y1, double z1, double x2, double y2, double z2, 
        double x3, double y3, double z3, double xe, double ye, double ze)
{
    MatrixXd matResult;

    matResult = ListMat(ListVec((1.*pow(x2 - xe,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (4.*(x1 - x3)*(x2 - xe)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (4.*pow(x1 - x3,2)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (1.*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(x2 - xe)*(y2 - ye))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(y1 - y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x1 - x3)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (4.*(x1 - x3)*(y1 - y3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(x2 - xe)*(z2 - ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x1 - x3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (4.*(x1 - x3)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(x1 - x3)*(x2 - xe))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*pow(x2 - xe,2)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*pow(x1 - x3,2)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (1.*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x1 - x3)*(x2 - xe)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(x2 - xe)*(y1 - y3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(x1 - x3)*(y1 - y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x1 - x3)*(y2 - ye)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(x2 - xe)*(z1 - z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x1 - x3)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(x1 - x3)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(x2 - xe)*(-x2 + xe))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x1 - x3)*(x2 - xe)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x1 - x3)*(-x2 + xe)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (4.*pow(x1 - x3,2)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (1.*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(x2 - xe)*(-y2 + ye))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(y1 - y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x1 - x3)*(-y2 + ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (4.*(x1 - x3)*(y1 - y3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(x2 - xe)*(-z2 + ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x1 - x3)*(-z2 + ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (4.*(x1 - x3)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(-x1 + x3)*(x2 - xe))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*pow(x2 - xe,2)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(x1 - x3)*(-x1 + x3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (1.*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x1 - x3)*(x2 - xe)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(x2 - xe)*(-y1 + y3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(x1 - x3)*(-y1 + y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x1 - x3)*(y2 - ye)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(x2 - xe)*(-z1 + z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x1 - x3)*(-z1 + z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(x1 - x3)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2))),
   ListVec((1.*(x2 - xe)*(y2 - ye))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(y1 - y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x1 - x3)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (4.*(x1 - x3)*(y1 - y3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*pow(y2 - ye,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (4.*(y1 - y3)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (4.*pow(y1 - y3,2)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (1.*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(y2 - ye)*(z2 - ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y2 - ye)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y1 - y3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (4.*(y1 - y3)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(x1 - x3)*(y2 - ye))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(x1 - x3)*(y1 - y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(y1 - y3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(y1 - y3)*(y2 - ye))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*pow(y2 - ye,2)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*pow(y1 - y3,2)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (1.*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y1 - y3)*(y2 - ye)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(y2 - ye)*(z1 - z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y1 - y3)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y2 - ye)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(y1 - y3)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-x2 + xe)*(y2 - ye))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(-x2 + xe)*(y1 - y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x1 - x3)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (4.*(x1 - x3)*(y1 - y3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(y2 - ye)*(-y2 + ye))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y1 - y3)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y1 - y3)*(-y2 + ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (4.*pow(y1 - y3,2)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (1.*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(y2 - ye)*(-z2 + ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y2 - ye)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y1 - y3)*(-z2 + ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (4.*(y1 - y3)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(-x1 + x3)*(y2 - ye))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(-x1 + x3)*(y1 - y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(y1 - y3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-y1 + y3)*(y2 - ye))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*pow(y2 - ye,2)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(y1 - y3)*(-y1 + y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (1.*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y1 - y3)*(y2 - ye)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(y2 - ye)*(-z1 + z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y1 - y3)*(-z1 + z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y2 - ye)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(y1 - y3)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2))),
   ListVec((1.*(x2 - xe)*(z2 - ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x1 - x3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (4.*(x1 - x3)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(y2 - ye)*(z2 - ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y2 - ye)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y1 - y3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (4.*(y1 - y3)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*pow(z2 - ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (4.*(z1 - z3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (1.*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (4.*pow(z1 - z3,2)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(x1 - x3)*(z2 - ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x1 - x3)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(x2 - xe)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(y1 - y3)*(z2 - ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y1 - y3)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y2 - ye)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(y2 - ye)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(z1 - z3)*(z2 - ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (1.*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*pow(z1 - z3,2)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*pow(z2 - ze,2)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(z1 - z3)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-x2 + xe)*(z2 - ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(-x2 + xe)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x1 - x3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (4.*(x1 - x3)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(-y2 + ye)*(z2 - ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(-y2 + ye)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y1 - y3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (4.*(y1 - y3)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(z2 - ze)*(-z2 + ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(z1 - z3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(z1 - z3)*(-z2 + ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (1.*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (4.*pow(z1 - z3,2)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(-x1 + x3)*(z2 - ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(-x1 + x3)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(x2 - xe)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-y1 + y3)*(z2 - ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(-y1 + y3)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y2 - ye)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(y2 - ye)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-z1 + z3)*(z2 - ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (1.*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(z1 - z3)*(-z1 + z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*pow(z2 - ze,2)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(z1 - z3)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2))),
   ListVec((1.*(x1 - x3)*(x2 - xe))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*pow(x2 - xe,2)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*pow(x1 - x3,2)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (1.*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x1 - x3)*(x2 - xe)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(x1 - x3)*(y2 - ye))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(x1 - x3)*(y1 - y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(y1 - y3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(x1 - x3)*(z2 - ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x1 - x3)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(x2 - xe)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*pow(x1 - x3,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (4.*(x1 - x3)*(x2 - xe)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (4.*pow(x2 - xe,2)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)) - 
     (1.*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(x1 - x3)*(y1 - y3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(y1 - y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(x1 - x3)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (4.*(x2 - xe)*(y2 - ye)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)),
    (1.*(x1 - x3)*(z1 - z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(x1 - x3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (4.*(x2 - xe)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)),
    (1.*(x1 - x3)*(-x2 + xe))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(-x2 + xe)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*pow(x1 - x3,2)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (1.*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x1 - x3)*(x2 - xe)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(x1 - x3)*(-y2 + ye))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(-y2 + ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(x1 - x3)*(y1 - y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(y1 - y3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(x1 - x3)*(-z2 + ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x1 - x3)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(-z2 + ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(x2 - xe)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(x1 - x3)*(-x1 + x3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x1 - x3)*(x2 - xe)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(-x1 + x3)*(x2 - xe)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (4.*pow(x2 - xe,2)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)) + 
     (1.*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(x1 - x3)*(-y1 + y3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(-y1 + y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(x1 - x3)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (4.*(x2 - xe)*(y2 - ye)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)),
    (1.*(x1 - x3)*(-z1 + z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(-z1 + z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(x1 - x3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (4.*(x2 - xe)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3))),
   ListVec((1.*(x2 - xe)*(y1 - y3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(x1 - x3)*(y1 - y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x1 - x3)*(y2 - ye)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(y1 - y3)*(y2 - ye))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*pow(y2 - ye,2)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*pow(y1 - y3,2)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (1.*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y1 - y3)*(y2 - ye)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(y1 - y3)*(z2 - ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y1 - y3)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y2 - ye)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(y2 - ye)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(x1 - x3)*(y1 - y3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(y1 - y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(x1 - x3)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (4.*(x2 - xe)*(y2 - ye)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)),
    (1.*pow(y1 - y3,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (4.*(y1 - y3)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (4.*pow(y2 - ye,2)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)) - 
     (1.*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(y1 - y3)*(z1 - z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y2 - ye)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(y1 - y3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (4.*(y2 - ye)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)),
    (1.*(-x2 + xe)*(y1 - y3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(-x2 + xe)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(x1 - x3)*(y1 - y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x1 - x3)*(y2 - ye)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(y1 - y3)*(-y2 + ye))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y2 - ye)*(-y2 + ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*pow(y1 - y3,2)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (1.*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y1 - y3)*(y2 - ye)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(y1 - y3)*(-z2 + ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y1 - y3)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y2 - ye)*(-z2 + ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(y2 - ye)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-x1 + x3)*(y1 - y3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(y1 - y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(-x1 + x3)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (4.*(x2 - xe)*(y2 - ye)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)),
    (1.*(y1 - y3)*(-y1 + y3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y1 - y3)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(-y1 + y3)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (4.*pow(y2 - ye,2)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)) + 
     (1.*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(y1 - y3)*(-z1 + z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y2 - ye)*(-z1 + z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(y1 - y3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (4.*(y2 - ye)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3))),
   ListVec((1.*(x2 - xe)*(z1 - z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x1 - x3)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(x1 - x3)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(y2 - ye)*(z1 - z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y1 - y3)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y2 - ye)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(y1 - y3)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(z1 - z3)*(z2 - ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (1.*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*pow(z1 - z3,2)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*pow(z2 - ze,2)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(z1 - z3)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(x1 - x3)*(z1 - z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(x1 - x3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (4.*(x2 - xe)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)),
    (1.*(y1 - y3)*(z1 - z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y2 - ye)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(y1 - y3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (4.*(y2 - ye)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)),
    (1.*pow(z1 - z3,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (4.*(z1 - z3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (1.*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (4.*pow(z2 - ze,2)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)),
    (1.*(-x2 + xe)*(z1 - z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x1 - x3)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(-x2 + xe)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(x1 - x3)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-y2 + ye)*(z1 - z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y1 - y3)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(-y2 + ye)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(y1 - y3)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(z1 - z3)*(-z2 + ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (1.*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*pow(z1 - z3,2)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(z2 - ze)*(-z2 + ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(z1 - z3)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-x1 + x3)*(z1 - z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(-x1 + x3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (4.*(x2 - xe)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)),
    (1.*(-y1 + y3)*(z1 - z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y2 - ye)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(-y1 + y3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (4.*(y2 - ye)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)),
    (1.*(z1 - z3)*(-z1 + z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(z1 - z3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(-z1 + z3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (1.*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (4.*pow(z2 - ze,2)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3))),
   ListVec((1.*(x2 - xe)*(-x2 + xe))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x1 - x3)*(x2 - xe)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x1 - x3)*(-x2 + xe)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (4.*pow(x1 - x3,2)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (1.*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(-x2 + xe)*(y2 - ye))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(-x2 + xe)*(y1 - y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x1 - x3)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (4.*(x1 - x3)*(y1 - y3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(-x2 + xe)*(z2 - ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(-x2 + xe)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x1 - x3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (4.*(x1 - x3)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(x1 - x3)*(-x2 + xe))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(-x2 + xe)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*pow(x1 - x3,2)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (1.*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x1 - x3)*(x2 - xe)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-x2 + xe)*(y1 - y3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(-x2 + xe)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(x1 - x3)*(y1 - y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x1 - x3)*(y2 - ye)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-x2 + xe)*(z1 - z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x1 - x3)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(-x2 + xe)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(x1 - x3)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*pow(-x2 + xe,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (4.*(x1 - x3)*(-x2 + xe)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (4.*pow(x1 - x3,2)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (1.*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(-x2 + xe)*(-y2 + ye))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(-x2 + xe)*(y1 - y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x1 - x3)*(-y2 + ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (4.*(x1 - x3)*(y1 - y3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(-x2 + xe)*(-z2 + ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(-x2 + xe)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x1 - x3)*(-z2 + ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (4.*(x1 - x3)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(-x1 + x3)*(-x2 + xe))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(-x2 + xe)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(x1 - x3)*(-x1 + x3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (1.*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x1 - x3)*(x2 - xe)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-x2 + xe)*(-y1 + y3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(-x2 + xe)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(x1 - x3)*(-y1 + y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x1 - x3)*(y2 - ye)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-x2 + xe)*(-z1 + z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x1 - x3)*(-z1 + z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(-x2 + xe)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(x1 - x3)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2))),
   ListVec((1.*(x2 - xe)*(-y2 + ye))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(y1 - y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x1 - x3)*(-y2 + ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (4.*(x1 - x3)*(y1 - y3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(y2 - ye)*(-y2 + ye))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y1 - y3)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y1 - y3)*(-y2 + ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (4.*pow(y1 - y3,2)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (1.*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(-y2 + ye)*(z2 - ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(-y2 + ye)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y1 - y3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (4.*(y1 - y3)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(x1 - x3)*(-y2 + ye))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(-y2 + ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(x1 - x3)*(y1 - y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(y1 - y3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(y1 - y3)*(-y2 + ye))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y2 - ye)*(-y2 + ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*pow(y1 - y3,2)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (1.*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y1 - y3)*(y2 - ye)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-y2 + ye)*(z1 - z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y1 - y3)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(-y2 + ye)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(y1 - y3)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-x2 + xe)*(-y2 + ye))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(-x2 + xe)*(y1 - y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x1 - x3)*(-y2 + ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (4.*(x1 - x3)*(y1 - y3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*pow(-y2 + ye,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (4.*(y1 - y3)*(-y2 + ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (4.*pow(y1 - y3,2)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (1.*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(-y2 + ye)*(-z2 + ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(-y2 + ye)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y1 - y3)*(-z2 + ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (4.*(y1 - y3)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(-x1 + x3)*(-y2 + ye))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(-y2 + ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(-x1 + x3)*(y1 - y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(y1 - y3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-y1 + y3)*(-y2 + ye))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y2 - ye)*(-y2 + ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(y1 - y3)*(-y1 + y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (1.*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y1 - y3)*(y2 - ye)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-y2 + ye)*(-z1 + z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y1 - y3)*(-z1 + z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(-y2 + ye)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(y1 - y3)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2))),
   ListVec((1.*(x2 - xe)*(-z2 + ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x1 - x3)*(-z2 + ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (4.*(x1 - x3)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(y2 - ye)*(-z2 + ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y2 - ye)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y1 - y3)*(-z2 + ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (4.*(y1 - y3)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(z2 - ze)*(-z2 + ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(z1 - z3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(z1 - z3)*(-z2 + ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (1.*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (4.*pow(z1 - z3,2)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(x1 - x3)*(-z2 + ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x1 - x3)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(-z2 + ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(x2 - xe)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(y1 - y3)*(-z2 + ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y1 - y3)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y2 - ye)*(-z2 + ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(y2 - ye)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(z1 - z3)*(-z2 + ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (1.*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*pow(z1 - z3,2)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(z2 - ze)*(-z2 + ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(z1 - z3)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-x2 + xe)*(-z2 + ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(-x2 + xe)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x1 - x3)*(-z2 + ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (4.*(x1 - x3)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(-y2 + ye)*(-z2 + ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(-y2 + ye)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y1 - y3)*(-z2 + ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (4.*(y1 - y3)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*pow(-z2 + ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (4.*(z1 - z3)*(-z2 + ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (1.*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (4.*pow(z1 - z3,2)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),3)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))),
    (1.*(-x1 + x3)*(-z2 + ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(-x1 + x3)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(-z2 + ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(x2 - xe)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-y1 + y3)*(-z2 + ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(-y1 + y3)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y2 - ye)*(-z2 + ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(y2 - ye)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-z1 + z3)*(-z2 + ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (1.*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(z1 - z3)*(-z1 + z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(z2 - ze)*(-z2 + ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(z1 - z3)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2))),
   ListVec((1.*(-x1 + x3)*(x2 - xe))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*pow(x2 - xe,2)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(x1 - x3)*(-x1 + x3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (1.*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x1 - x3)*(x2 - xe)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-x1 + x3)*(y2 - ye))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(-x1 + x3)*(y1 - y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(y1 - y3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-x1 + x3)*(z2 - ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(-x1 + x3)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(x2 - xe)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(x1 - x3)*(-x1 + x3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x1 - x3)*(x2 - xe)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(-x1 + x3)*(x2 - xe)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (4.*pow(x2 - xe,2)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)) + 
     (1.*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-x1 + x3)*(y1 - y3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(y1 - y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(-x1 + x3)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (4.*(x2 - xe)*(y2 - ye)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)),
    (1.*(-x1 + x3)*(z1 - z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(-x1 + x3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (4.*(x2 - xe)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)),
    (1.*(-x1 + x3)*(-x2 + xe))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(-x2 + xe)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(x1 - x3)*(-x1 + x3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (1.*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x1 - x3)*(x2 - xe)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-x1 + x3)*(-y2 + ye))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(-y2 + ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(-x1 + x3)*(y1 - y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(y1 - y3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-x1 + x3)*(-z2 + ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(-x1 + x3)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(-z2 + ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(x2 - xe)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*pow(-x1 + x3,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (4.*(-x1 + x3)*(x2 - xe)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (4.*pow(x2 - xe,2)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)) - 
     (1.*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-x1 + x3)*(-y1 + y3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(-y1 + y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(-x1 + x3)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (4.*(x2 - xe)*(y2 - ye)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)),
    (1.*(-x1 + x3)*(-z1 + z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(-z1 + z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(-x1 + x3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (4.*(x2 - xe)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3))),
   ListVec((1.*(x2 - xe)*(-y1 + y3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(x1 - x3)*(-y1 + y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x1 - x3)*(y2 - ye)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-y1 + y3)*(y2 - ye))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*pow(y2 - ye,2)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(y1 - y3)*(-y1 + y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (1.*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y1 - y3)*(y2 - ye)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-y1 + y3)*(z2 - ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(-y1 + y3)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y2 - ye)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(y2 - ye)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(x1 - x3)*(-y1 + y3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(-y1 + y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(x1 - x3)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (4.*(x2 - xe)*(y2 - ye)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)),
    (1.*(y1 - y3)*(-y1 + y3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y1 - y3)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(-y1 + y3)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (4.*pow(y2 - ye,2)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)) + 
     (1.*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-y1 + y3)*(z1 - z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y2 - ye)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(-y1 + y3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (4.*(y2 - ye)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)),
    (1.*(-x2 + xe)*(-y1 + y3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(-x2 + xe)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(x1 - x3)*(-y1 + y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x1 - x3)*(y2 - ye)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-y1 + y3)*(-y2 + ye))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y2 - ye)*(-y2 + ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(y1 - y3)*(-y1 + y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (1.*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y1 - y3)*(y2 - ye)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-y1 + y3)*(-z2 + ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(-y1 + y3)*(z1 - z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y2 - ye)*(-z2 + ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(y2 - ye)*(z1 - z3)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-x1 + x3)*(-y1 + y3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(-y1 + y3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(-x1 + x3)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (4.*(x2 - xe)*(y2 - ye)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)),
    (1.*pow(-y1 + y3,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (4.*(-y1 + y3)*(y2 - ye)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (4.*pow(y2 - ye,2)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)) - 
     (1.*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-y1 + y3)*(-z1 + z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y2 - ye)*(-z1 + z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(-y1 + y3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (4.*(y2 - ye)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3))),
   ListVec((1.*(x2 - xe)*(-z1 + z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x1 - x3)*(-z1 + z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(x1 - x3)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(y2 - ye)*(-z1 + z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y1 - y3)*(-z1 + z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y2 - ye)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(y1 - y3)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-z1 + z3)*(z2 - ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (1.*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(z1 - z3)*(-z1 + z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*pow(z2 - ze,2)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(z1 - z3)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(x1 - x3)*(-z1 + z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(x2 - xe)*(-z1 + z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(x1 - x3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (4.*(x2 - xe)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)),
    (1.*(y1 - y3)*(-z1 + z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) - 
     (2.*(y2 - ye)*(-z1 + z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(y1 - y3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (4.*(y2 - ye)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)),
    (1.*(z1 - z3)*(-z1 + z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(z1 - z3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (2.*(-z1 + z3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (1.*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (4.*pow(z2 - ze,2)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)),
    (1.*(-x2 + xe)*(-z1 + z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x1 - x3)*(-z1 + z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(-x2 + xe)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(x1 - x3)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-y2 + ye)*(-z1 + z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y1 - y3)*(-z1 + z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(-y2 + ye)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(y1 - y3)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-z1 + z3)*(-z2 + ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (1.*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(z1 - z3)*(-z1 + z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(z2 - ze)*(-z2 + ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(z1 - z3)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      (pow(pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2),2)*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)),
    (1.*(-x1 + x3)*(-z1 + z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(x2 - xe)*(-z1 + z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(-x1 + x3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (4.*(x2 - xe)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)),
    (1.*(-y1 + y3)*(-z1 + z3))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (2.*(y2 - ye)*(-z1 + z3)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (2.*(-y1 + y3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (4.*(y2 - ye)*(z2 - ze)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3)),
    (1.*pow(-z1 + z3,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        (pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2))) + 
     (4.*(-z1 + z3)*(z2 - ze)*(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + 
          y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) - 
     (1.*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - y1*ye + y3*ye + z1*z2 - 
          z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),2)) + 
     (4.*pow(z2 - ze,2)*pow(-(x2*x3) + x1*(x2 - xe) + x3*xe + y1*y2 - y2*y3 - 
          y1*ye + y3*ye + z1*z2 - z2*z3 - z1*ze + z3*ze,2))/
      ((pow(x1 - x3,2) + pow(y1 - y3,2) + pow(z1 - z3,2))*
        pow(pow(x2 - xe,2) + pow(y2 - ye,2) + pow(z2 - ze,2),3))));

    return matResult;
}

VectorXd elasticAngleBound::ListVec(double a1, double a2, double a3, 
    double a4, double a5, double a6, double a7, double a8, double a9, 
    double a10, double a11, double a12)
{
    VectorXd vecResult;

    vecResult.setZero(12, 1);

    vecResult(0)  = a1;
    vecResult(1)  = a2;
    vecResult(2)  = a3;
    vecResult(3)  = a4;
    vecResult(4)  = a5;
    vecResult(5)  = a6;
    vecResult(6)  = a7;
    vecResult(7)  = a8;
    vecResult(8)  = a9;
    vecResult(9)  = a10;
    vecResult(10) = a11;
    vecResult(11) = a12;

    return vecResult;
}

MatrixXd elasticAngleBound::ListMat(VectorXd a1, VectorXd a2, VectorXd a3, 
    VectorXd a4, VectorXd a5, VectorXd a6, VectorXd a7, VectorXd a8,
    VectorXd a9, VectorXd a10, VectorXd a11, VectorXd a12)
{
    MatrixXd matResult;

    matResult.setZero(12, 12);

    matResult.col(0)  = a1;
    matResult.col(1)  = a2;
    matResult.col(2)  = a3;
    matResult.col(3)  = a4;
    matResult.col(4)  = a5;
    matResult.col(5)  = a6;
    matResult.col(6)  = a7;
    matResult.col(7)  = a8;
    matResult.col(8)  = a9;
    matResult.col(9)  = a10;
    matResult.col(10) = a11;
    matResult.col(11) = a12;

    return matResult;
}