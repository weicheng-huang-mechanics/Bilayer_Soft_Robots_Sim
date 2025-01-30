#include "externalGravityForce.h"
#include "elasticPlate.h"
#include "timeStepper.h"

externalGravityForce::externalGravityForce(const std::shared_ptr<elasticPlate> m_plate,
									       const std::shared_ptr<timeStepper> m_stepper,
										   Vector3d m_gVector)
										   : BaseForce(m_plate, m_stepper), gVector(m_gVector)
{
	setGravity();
}

externalGravityForce::~externalGravityForce() = default;

void externalGravityForce::computeForceAndJacobian()
{
	computeFg();
	computeJg();
}


void externalGravityForce::computeFg()
{
	for (int i=0; i < plate->ndof; i++)
	{
		stepper->addForce(i, -massGravity[i]); // subtracting gravity force
	}	
}

void externalGravityForce::computeJg()
{
	;
}

void externalGravityForce::setGravity()
{
	massGravity = VectorXd::Zero(plate->ndof);
	
	for (int i = 0; i < plate->nv; i++)
	{
		for (int k = 0; k < 3; k++)
		{
			int ind = 3 * i + k;
			massGravity[ind] = gVector[k] * plate->massArray[ind];
		}
	}
}

void externalGravityForce::setFirstJacobian()
{
	for (int i = 0; i < plate->nv; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			int ind = 3 * i + j;
			stepper->addJacobian(ind, ind, 1);
		}
	}
}