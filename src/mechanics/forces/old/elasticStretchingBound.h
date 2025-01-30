#ifndef ELASTICSTRETCHINGBOUND_H
#define ELASTICSTRETCHINGBOUND_H

#include "eigenIncludes.h"
#include "elasticPlate.h"
#include "timeStepper.h"

class elasticStretchingBound
{
public:
	elasticStretchingBound(elasticPlate &m_plate, timeStepper &m_stepper);
	~elasticStretchingBound();
	void computeFs();
	void computeJs();
    
    void setFirstJacobian();
    
private:
	elasticPlate *plate;
    timeStepper *stepper;

    VectorXi localDOF;

    VectorXd flocal;
    MatrixXd Jss;

    double EA;

    Matrix3d Id3;

    void localForce();
    void localJacobian();
};

#endif
