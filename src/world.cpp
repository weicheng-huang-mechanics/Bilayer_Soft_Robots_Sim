#include "world.h"

world::world()
{
	;
}

world::world(setInput &m_inputData)
{
	render = m_inputData.GetBoolOpt("render");				
	saveData = m_inputData.GetBoolOpt("saveData");			
	deltaTime = m_inputData.GetScalarOpt("deltaTime");     
	totalTime = m_inputData.GetScalarOpt("totalTime");    
	YoungM = m_inputData.GetScalarOpt("YoungM");
	density = m_inputData.GetScalarOpt("density");
	rodRadius = m_inputData.GetScalarOpt("rodRadius");
	Possion = m_inputData.GetScalarOpt("Possion");
	stol = m_inputData.GetScalarOpt("stol");
	forceTol = m_inputData.GetScalarOpt("forceTol");
	scaleRendering = m_inputData.GetScalarOpt("scaleRendering");
	maxIter = m_inputData.GetIntOpt("maxIter");
	gVector = m_inputData.GetVecOpt("gVector");
	viscosity = m_inputData.GetScalarOpt("viscosity");
}

world::~world()
{
	;
}

bool world::isRender()
{
	return render;
}

void world::OpenFile(ofstream &outfile)
{
	if (saveData==false) return;
	
	int systemRet = system("mkdir datafiles"); //make the directory
	if(systemRet == -1)
	{
		cout << "Error in creating directory\n";
	}

	// Open an input file named after the current time
	ostringstream name;
	name.precision(6);
	name << fixed;

    name << "datafiles/simDER";
    name << ".txt";

	outfile.open(name.str().c_str());
	outfile.precision(10);	
}

void world::CloseFile(ofstream &outfile)
{
	if (saveData==false) 
	{
		return;
	}
}

void world::CoutData(ofstream &outfile)
{
	if (saveData==false) 
	{
		return;
	}

	if ( timeStep % 10 != 0)
	{
		return;
	}

	if (timeStep == Nstep)
	{
		for (int i = 0; i < plate->nv; i++)
		{
			Vector3d xCurrent = plate->getVertex(i);

			outfile << xCurrent(0) << " " << xCurrent(1) << " " << xCurrent(2) << endl;
		}
	}
}

void world::setPlateStepper()
{
	// Create the plate 
	plate = std::make_shared<elasticPlate>(YoungM, density, rodRadius, Possion, deltaTime);
	plateBoundaryCondition();
	plate->setup();
	stepper = std::make_shared<timeStepper>(*plate);

	forces.clear();
	// set up force
	m_inertialForce = std::make_shared<inertialForce>(plate, stepper);
	m_gravityForce = std::make_shared<externalGravityForce>(plate, stepper, gVector);
	m_dampingForce = std::make_shared<dampingForce>(plate, stepper, viscosity);
	m_stretchForce = std::make_shared<elasticStretchingForce>(plate, stepper);
	m_bendingForce = std::make_shared<elasticBendingForce>(plate, stepper);
	m_twistingForce = std::make_shared<elasticTwistingForce>(plate, stepper);

	// set up constraints
	m_elasticBendingBound = std::make_shared<elasticBendingBound>(plate, stepper);
	m_elasticTwistingBound = std::make_shared<elasticTwistingBound>(plate, stepper);
	m_constrainedForce = std::make_shared<constrainedForce>(plate, stepper, 10000 * plate->EA, 0.005);	

	forces.push_back(m_inertialForce);
	forces.push_back(m_gravityForce);
	forces.push_back(m_dampingForce);
	forces.push_back(m_stretchForce);
	forces.push_back(m_bendingForce);
	forces.push_back(m_twistingForce);
	forces.push_back(m_elasticBendingBound);
	forces.push_back(m_elasticTwistingBound);
	forces.push_back(m_constrainedForce);
	
	
	plate->updateTimeStep();

	for (auto &force : forces)
	{
		force->setFirstJacobian();
	}
	stepper->first_time_PARDISO_setup();


	// time step 
	Nstep = totalTime / deltaTime;
	timeStep = 0;
	currentTime = 0.0;
}

void world::plateBoundaryCondition()
{
	plate->setVertexBoundaryCondition(plate->getVertex(0), 0);
	plate->setVertexBoundaryCondition(plate->getVertex(1), 1);
	plate->setThetaBoundaryCondition(plate->getTheta(0), 0);

	for (int i = 0; i < plate->v_edgeElement.size(); i++)
	{
		//plate->setThetaBoundaryCondition(plate->getTheta(i), i);
	}
}

void world::updateBoundaryCondition()
{
	if (currentTime < 10.0)
	{
		for (int i = 0; i < 196/2; i++)
		{
			plate->v_bendingElement[i].unTwist = currentTime * 0.01;
			plate->v_bendingElement[i].kappaBar(0) = currentTime * 0.01 * sin(6 * M_PI * i / 98);
			plate->v_bendingElement[i].kappaBar(1) = currentTime * 0.01 * cos(6 * M_PI * i / 98);
		}
	}
}

void world::updateTimeStep()
{
	plate->updateGuess(); // x = x0 + u * dt

	updateBoundaryCondition();
	updateEachStep();
	plate->updateTimeStep();


	if (render) 
	{
		cout << "time: " << currentTime << " ";
	}

	currentTime += deltaTime;
		
	timeStep++;
}


void world::updateEachStep()
{
	double normf = forceTol * 10.0;
	double normf0 = 0;
	
	bool solved = false;
	
	int iter = 0;
		
	while (solved == false)
	{
		plate->prepareForIteration();

		stepper->setZero();

		for (auto &force : forces)
		{
			force->computeForceAndJacobian();
		}


		normf = stepper->GlobalForceVec.norm();

		if (iter == 0) 
		{
			normf0 = normf;
		}
		
		if (normf <= forceTol)
		{
			solved = true;
		}
		else if(iter > 0 && normf <= normf0 * stol)
		{
			solved = true;
		}

		normf = 0.0;
		
		if (solved == false)
		{
			stepper->integrator(); // Solve equations of motion
			plate->updateNewtonMethod(stepper->GlobalMotionVec); // new q = old q + Delta q
			iter++;
		}

		if (iter > maxIter)
		{
			cout << "Error. Could not converge. Exiting.\n";
			break;
		}
	}

	if (render)
	{
		cout << "iter " << iter << endl;
	}
}

int world::simulationRunning()
{
	if (timeStep < Nstep) 
	{
		return 1;
	}
	else 
	{
		return -1;
	}
}

Vector3d world::getScaledCoordinate(int i, int j)
{
	Vector3d xCurrent;
	
	if (j == 0)
	{
		xCurrent = plate->v_edgeElement[i].x_1 * scaleRendering;
	}
	if (j == 1)
	{
		xCurrent = plate->v_edgeElement[i].x_2 * scaleRendering;
	}

	return xCurrent;
}

int world::numStretchingPair()
{
	return plate->edgeNum;
}

int world::numBendingPair()
{
	return plate->v_bendingElement.size();
}

Vector3d world::getBoundaryCoordination_left(int i, int j)
{
	Vector3d xCurrent;

	if (j == 0)
	{
		xCurrent = (plate->v_bendingElement[i].x_1 - 0.01 * plate->v_bendingElement[i].m_12) * scaleRendering;
	}
	if (j == 1)
	{
		xCurrent = (plate->v_bendingElement[i].x_3 - 0.01 * plate->v_bendingElement[i].m_22) * scaleRendering;
	}

	return xCurrent;
	
}

Vector3d world::getBoundaryCoordination_right(int i, int j)
{
	Vector3d xCurrent;

	if (j == 0)
	{
		xCurrent = (plate->v_bendingElement[i].x_1 + 0.01 * plate->v_bendingElement[i].m_12) * scaleRendering;
	}
	if (j == 1)
	{
		xCurrent = (plate->v_bendingElement[i].x_3 + 0.01 * plate->v_bendingElement[i].m_22) * scaleRendering;
	}

	return xCurrent;
}
