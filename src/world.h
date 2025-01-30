#ifndef WORLD_H
#define WORLD_H

#include "eigenIncludes.h"

// include input file and option
#include "setInput.h"

// include elastic Plate class
#include "elasticPlate.h"

// include time stepper
#include "timeStepper.h"

// include force
#include "mechanics/forces/inertialForce.h"
#include "mechanics/forces/externalGravityForce.h"
#include "mechanics/forces/dampingForce.h"
#include "mechanics/forces/elasticStretchingForce.h"
#include "mechanics/forces/elasticBendingForce.h"
#include "mechanics/forces/elasticTwistingForce.h"

// include constraints
#include "mechanics/constraints/elasticBendingBound.h"
#include "mechanics/constraints/elasticTwistingBound.h"
#include "mechanics/constraints/constrainedForce.h"


class world
{
public:
	world();
	world(setInput &m_inputData);
	~world();
	
	bool isRender();
	
	// file output
	void OpenFile(ofstream &outfile);
	void CloseFile(ofstream &outfile);
	void CoutData(ofstream &outfile);

	void setPlateStepper();

	void updateTimeStep();

	int simulationRunning();

	int numStretchingPair();
	Vector3d getScaledCoordinate(int i, int j);

	int numBendingPair();
	Vector3d getBoundaryCoordination_left(int i, int j);
	Vector3d getBoundaryCoordination_right(int i, int j);

private:

	// physical parameters
	bool render;
	bool saveData;
	double deltaTime;
	double totalTime;
	double YoungM;
	double density;
	double Possion;
	double stol;
	double forceTol;
	double scaleRendering;
	int maxIter;
	Vector3d gVector;
	double viscosity;
	double rodRadius;
	
	int Nstep;
	int timeStep;

	double characteristicForce;

	double currentTime;

	void plateBoundaryCondition();
	void updateBoundaryCondition();

	shared_ptr<elasticPlate> plate;
	shared_ptr<timeStepper> stepper;

	// force
	shared_ptr<inertialForce> m_inertialForce;
	shared_ptr<externalGravityForce> m_gravityForce;
	shared_ptr<dampingForce> m_dampingForce;
	shared_ptr<elasticStretchingForce> m_stretchForce;
	shared_ptr<elasticBendingForce> m_bendingForce;
	shared_ptr<elasticTwistingForce> m_twistingForce;

	shared_ptr<elasticBendingBound> m_elasticBendingBound;	
	shared_ptr<elasticTwistingBound> m_elasticTwistingBound;
	shared_ptr<constrainedForce> m_constrainedForce;

	vector<shared_ptr<BaseForce>> forces;

	void updateEachStep();
};

#endif
