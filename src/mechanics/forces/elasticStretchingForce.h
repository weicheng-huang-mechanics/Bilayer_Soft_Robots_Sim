#ifndef ELASTICSTRETCHINGFORCE_H
#define ELASTICSTRETCHINGFORCE_H

#include "mechanics/base_force.h"

class elasticStretchingForce : public BaseForce
{
public:
	elasticStretchingForce(const std::shared_ptr<elasticPlate> m_plate, const std::shared_ptr<timeStepper> m_stepper); 
	~elasticStretchingForce() override;

    void computeForceAndJacobian() override;
    void setFirstJacobian() override;
    
private:
	void computeFs();
	void computeJs();

    VectorXi localDOF;

    VectorXd flocal;
    MatrixXd Jss;

    double EA;

    Matrix3d Id3;

    void localForce();
    void localJacobian();
};

#endif
