#ifndef ELASTICANGLEBOUND_H
#define ELASTICANGLEBOUND_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class elasticAngleBound
{
public:
	elasticAngleBound(elasticPlate &m_plate, timeStepper &m_stepper);
	~elasticAngleBound();
	void computeFa();
	void computeJa();
    
    void setFirstJacobian();
    
private:
	elasticPlate *plate;
    timeStepper *stepper;

    double EA;

    void localForce();
    void localJacobian();

    MatrixXd ListMat(VectorXd a1, VectorXd a2, VectorXd a3, 
    VectorXd a4, VectorXd a5, VectorXd a6, VectorXd a7, VectorXd a8,
    VectorXd a9, VectorXd a10, VectorXd a11, VectorXd a12);
    VectorXd ListVec(double a1, double a2, double a3, 
    double a4, double a5, double a6, double a7, double a8, double a9, 
    double a10, double a11, double a12);

    VectorXd computeAngleForce(double x1, double y1, double z1, double x2, double y2, double z2, 
        double x3, double y3, double z3, double xe, double ye, double ze);
    MatrixXd computeAngleJacobian(double x1, double y1, double z1, double x2, double y2, double z2, 
        double x3, double y3, double z3, double xe, double ye, double ze);
};

#endif
