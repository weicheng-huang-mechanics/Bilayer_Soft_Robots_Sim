#ifndef INERTIALFORCE_H
#define INERTIALFORCE_H

#include "mechanics/base_force.h"

class inertialForce : public BaseForce
{
public:
	inertialForce(const std::shared_ptr<elasticPlate> m_plate, const std::shared_ptr<timeStepper> m_stepper);
	~inertialForce() override;

	void computeForceAndJacobian() override;
	void setFirstJacobian() override;

private:
	void computeFi();
	void computeJi();

	// support variables
    int ind1, ind2, mappedInd1, mappedInd2;	
    double f, jac;
};

#endif
