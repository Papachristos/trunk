 
/*************************************************************************
*  Copyright (C) 2014 by Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>     *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
	//TODO: Cubic Law on fractures (~~ l. 250 )
	//TODO: Add a " CrackWidth" factor (I guess, of the same order as the permeabFactor (??) 
	//TODO: Add a minimum value for the crackWidth --> Residual permeability (can be the one of the intact rock for instance) 


#ifdef FLOW_ENGINE

#include<yade/pkg/dem/JointedCohesiveFrictionalPM.hpp>

//keep this #ifdef for commited versions unless you really have stable version that should be compiled by default
//it will save compilation time for everyone else
//when you want it compiled, you can pass -DDFNFLOW to cmake, or just uncomment the following line
#define DFNFLOW
#ifdef DFNFLOW
#define TEMPLATE_FLOW_NAME DFNFlowEngineT
#include <yade/pkg/pfv/FlowEngine.hpp>

class DFNCellInfo : public FlowCellInfo
// Here we add every "extra attribute" we need for cells
{
	public:
	Real anotherVariable;
	bool crack;
// 	bool preExistingJoint;
	void anotherFunction() {};
	DFNCellInfo() : FlowCellInfo(),crack(false)  {}
};

// Here we add every extra attribute we need for the vertices
class DFNVertexInfo : public FlowVertexInfo {
	public:
	//same here if needed
};

typedef CGT::_Tesselation<CGT::TriangulationTypes<DFNVertexInfo,DFNCellInfo> > DFNTesselation;
// We add all the new/complementary things for FlowBoundingSphere.ipp in the DFNboundingSphere
// always declare public what we call from another class
#ifdef LINSOLV
class DFNBoundingSphere : public CGT::FlowBoundingSphereLinSolv<DFNTesselation>
#else
class DFNBoundingSphere : public CGT::FlowBoundingSphere<DFNTesselation>
#endif
{
public:
 // Added to the FlowBoundingSphere.ipp saveVtk to include fractured cells
 ///Old
//  void computePermeability()
// {
// 	if (debugOut)  cout << "----Computing_Permeability------" << endl;
// 	RTriangulation& Tri = T[currentTes].Triangulation();
// 	VSolidTot = 0, Vtotalissimo = 0, vPoral = 0, sSolidTot = 0, vTotalPorosity=0, vPoralPorosity=0;
// 	FiniteCellsIterator cellEnd = Tri.finite_cells_end();
// 
// 	CellHandle neighbourCell;
// 
// 	double k=0, distance = 0, radius = 0;
// 	int surfneg=0;
// 	int NEG=0, POS=0, pass=0;
// 
// 	bool ref = Tri.finite_cells_begin()->info().isvisited;
// 	Real meanK=0, STDEV=0, meanRadius=0, meanDistance=0;
// 	Real infiniteK=1e10;
// 
// 	for (VCellIterator cellIt=T[currentTes].cellHandles.begin(); cellIt!=T[currentTes].cellHandles.end(); cellIt++){
// 		CellHandle& cell = *cellIt;
// 		Point& p1 = cell->info();
// 		for (int j=0; j<4; j++) {
// 			neighbourCell = cell->neighbor(j);
// 			Point& p2 = neighbourCell->info();
// 			if (!Tri.is_infinite(neighbourCell) && (neighbourCell->info().isvisited==ref || computeAllCells)) {
// 				//compute and store the area of sphere-facet intersections for later use
// 				VertexHandle W [3];
// 				for (int kk=0; kk<3; kk++) {
// 					W[kk] = cell->vertex(facetVertices[j][kk]);
// 				}
// 				CGT::Sphere& v0 = W[0]->point();
// 				CGT::Sphere& v1 = W[1]->point();
// 				CGT::Sphere& v2 = W[2]->point();
// 
// 				cell->info().facetSphereCrossSections[j]=CVector(
// 				   W[0]->info().isFictious ? 0 : 0.5*v0.weight()*acos((v1-v0)*(v2-v0)/sqrt((v1-v0).squared_length()*(v2-v0).squared_length())),
// 				   W[1]->info().isFictious ? 0 : 0.5*v1.weight()*acos((v0-v1)*(v2-v1)/sqrt((v1-v0).squared_length()*(v2-v1).squared_length())),
// 				   W[2]->info().isFictious ? 0 : 0.5*v2.weight()*acos((v0-v2)*(v1-v2)/sqrt((v1-v2).squared_length()*(v2-v0).squared_length())));
// 
// 				pass+=1;
// 				CVector l = p1 - p2;
// 				distance = sqrt(l.squared_length());
// 				if (!rAverage) radius = 2* computeHydraulicRadius(cell, j);
// 				else radius = (computeEffectiveRadius(cell, j)+computeEquivalentRadius(cell,j))*0.5;
// 				if (radius<0) NEG++;
// 				else POS++;
// 				if (radius==0) {
// 					cout << "INS-INS PROBLEM!!!!!!!" << endl;
// 				}
// 				Real fluidArea=0;
// 				if (distance!=0) {
// 					if (minPermLength>0 && distanceCorrection) distance=max(minPermLength,distance);
// 					const CVector& Surfk = cell->info().facetSurfaces[j];
// 					Real area = sqrt(Surfk.squared_length());
// 					const CVector& crossSections = cell->info().facetSphereCrossSections[j];
// 					Real S0=0;
// 					S0=checkSphereFacetOverlap(v0,v1,v2);
// 					if (S0==0) S0=checkSphereFacetOverlap(v1,v2,v0);
// 					if (S0==0) S0=checkSphereFacetOverlap(v2,v0,v1);
// 					//take absolute value, since in rare cases the surface can be negative (overlaping spheres)
// 					fluidArea=abs(area-crossSections[0]-crossSections[1]-crossSections[2]+S0);
// 					cell->info().facetFluidSurfacesRatio[j]=fluidArea/area;
// 					k=(fluidArea * pow(radius,2)) / (8*viscosity*distance);
// 					 meanDistance += distance;
// 					 meanRadius += radius;
// 					 meanK +=  k*kFactor;
// 
// 				if (k<0 && debugOut) {surfneg+=1;
// 				cout<<"__ k<0 __"<<k<<" "<<" fluidArea "<<fluidArea<<" area "<<area<<" "<<crossSections[0]<<" "<<crossSections[1]<<" "<<crossSections[2] <<" "<<W[0]->info().id()<<" "<<W[1]->info().id()<<" "<<W[2]->info().id()<<" "<<p1<<" "<<p2<<" test "<<endl;}				     
// 				} else  {cout <<"infinite K1!"<<endl; k = infiniteK;}//Will be corrected in the next loop
// 
// 				(cell->info().kNorm())[j]= k*kFactor;
// 				if (!neighbourCell->info().isGhost) (neighbourCell->info().kNorm())[Tri.mirror_index(cell, j)]= (cell->info().kNorm())[j];
// 				
// 			}
// 		}
// 		cell->info().isvisited = !ref;
// 	}
// 	if (debugOut) cout<<"surfneg est "<<surfneg<<endl;
// 	meanK /= pass;
// 	meanRadius /= pass;
// 	meanDistance /= pass;
// 	Real globalK=kFactor*meanDistance*vPoral/(sSolidTot*8.*viscosity);//An approximate value of macroscopic permeability, for clamping local values below
// 	if (debugOut) {
// 		cout << "PassCompK = " << pass << endl;
// 		cout << "meanK = " << meanK << endl;
// 		cout << "globalK = " << globalK << endl;
// 		cout << "maxKdivKmean*globalK = " << maxKdivKmean*globalK << endl;
// 		cout << "minKdivKmean*globalK = " << minKdivKmean*globalK << endl;
// 		cout << "meanTubesRadius = " << meanRadius << endl;
// 		cout << "meanDistance = " << meanDistance << endl;
// 	}
// 	ref = Tri.finite_cells_begin()->info().isvisited;
// 	pass=0;
// 
// 	if (clampKValues) for (VCellIterator cellIt=T[currentTes].cellHandles.begin(); cellIt!=T[currentTes].cellHandles.end(); cellIt++){
// 		CellHandle& cell = *cellIt;
// 		for (int j=0; j<4; j++) {
// 			neighbourCell = cell->neighbor(j);
// 			if (!Tri.is_infinite(neighbourCell) && neighbourCell->info().isvisited==ref) {
// 				pass++;
// 				(cell->info().kNorm())[j] = max(minKdivKmean*globalK ,min((cell->info().kNorm())[j], maxKdivKmean*globalK));
// 				(neighbourCell->info().kNorm())[Tri.mirror_index(cell, j)]=(cell->info().kNorm())[j];
// 			}
// 		}
// 	}
// 	if (debugOut) cout << "PassKcorrect = " << pass << endl;
// 	if (debugOut) cout << "POS = " << POS << " NEG = " << NEG << " pass = " << pass << endl;
// 
// 	// A loop to compute the standard deviation of the local K distribution, and use it to include/exclude K values higher then (meanK +/- K_opt_factor*STDEV)
// 	if (meanKStat)
// 	{
// 		std::ofstream k_opt_file("k_stdev.txt" ,std::ios::out);
// 		ref = Tri.finite_cells_begin()->info().isvisited;
// 		pass=0;
// 		for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
// 			for (int j=0; j<4; j++) {
// 				neighbourCell = cell->neighbor(j);
// 				if (!Tri.is_infinite(neighbourCell) && neighbourCell->info().isvisited==ref) {
// 					pass++;
// 					STDEV += pow(((cell->info().kNorm())[j]-meanK),2);
// 				}
// 			}cell->info().isvisited = !ref;
// 		}
// 		STDEV = sqrt(STDEV/pass);
// 		if (debugOut) cout << "PassSTDEV = " << pass << endl << "STATISTIC K" << endl;
// 		double k_min = 0, k_max = meanK + KOptFactor*STDEV;
// 		cout << "Kmoy = " << meanK << " Standard Deviation = " << STDEV << endl<< "kmin = " << k_min << " kmax = " << k_max << endl;
// 		ref = Tri.finite_cells_begin()->info().isvisited;
// 		pass=0;
// 		for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
// 			for (int j=0; j<4; j++) {
// 				neighbourCell = cell->neighbor(j);
// 				if (!Tri.is_infinite(neighbourCell) && neighbourCell->info().isvisited==ref) {
// 					pass+=1;
// 					if ((cell->info().kNorm())[j]>k_max) {
// 						(cell->info().kNorm())[j]=k_max;
// 						(neighbourCell->info().kNorm())[Tri.mirror_index(cell, j)]= (cell->info().kNorm())[j];
// 					}
// 					k_opt_file << KOptFactor << " " << (cell->info().kNorm())[j] << endl;
// 				}
// 			}cell->info().isvisited=!ref;
// 		}
// 		if (debugOut) cout << "PassKopt = " << pass << endl;
// 	}
// 	if (debugOut) {
// 		FiniteVerticesIterator verticesEnd = Tri.finite_vertices_end();
// 		Real Vgrains = 0;
// 		int grains=0;
// 		for (FiniteVerticesIterator vIt = Tri.finite_vertices_begin(); vIt !=  verticesEnd; vIt++) {
// 			if (!vIt->info().isFictious && !vIt->info().isGhost) {
// 				grains +=1;
// 				Vgrains += 1.33333333 * M_PI * pow(vIt->point().weight(),1.5);}}
// 		cout<<grains<<"grains - " <<"vTotal = " << vTotal << " Vgrains = " << Vgrains << " vPoral1 = " << (vTotal-Vgrains) << endl;
// 		cout << "Vtotalissimo = " << Vtotalissimo/2 << " VSolidTot = " << VSolidTot/2 << " vPoral2 = " << vPoral/2  << " sSolidTot = " << sSolidTot << endl<< endl;
// 		if (!rAverage) cout << "------Hydraulic Radius is used for permeability computation------" << endl << endl;
// 		else cout << "------Average Radius is used for permeability computation------" << endl << endl;
// 		cout << "-----computed_Permeability-----" << endl;}
// }
/// Up to here

  void saveVtk(const char* folder)
  {
//     CGT::FlowBoundingSphere<DFNTesselation>::saveVtk(folder)
	RTriangulation& Tri = T[noCache?(!currentTes):currentTes].Triangulation();
        static unsigned int number = 0;
        char filename[80];
	mkdir(folder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        sprintf(filename,"%s/out_%d.vtk",folder,number++);
	int firstReal=-1;

	//count fictious vertices and cells
	vtkInfiniteVertices=vtkInfiniteCells=0;
 	FiniteCellsIterator cellEnd = Tri.finite_cells_end();
        for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != cellEnd; cell++) {
		bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal() && cell->vertex(2)->info().isReal()  && cell->vertex(3)->info().isReal();
		if (!isDrawable) vtkInfiniteCells+=1;
	}
	for (FiniteVerticesIterator v = Tri.finite_vertices_begin(); v != Tri.finite_vertices_end(); ++v) {
                if (!v->info().isReal()) vtkInfiniteVertices+=1;
                else if (firstReal==-1) firstReal=vtkInfiniteVertices;}

        basicVTKwritter vtkfile((unsigned int) Tri.number_of_vertices()-vtkInfiniteVertices, (unsigned int) Tri.number_of_finite_cells()-vtkInfiniteCells);

        vtkfile.open(filename,"test");

        vtkfile.begin_vertices();
        double x,y,z;
        for (FiniteVerticesIterator v = Tri.finite_vertices_begin(); v != Tri.finite_vertices_end(); ++v) {
		if (v->info().isReal()){
		x = (double)(v->point().point()[0]);
                y = (double)(v->point().point()[1]);
                z = (double)(v->point().point()[2]);
                vtkfile.write_point(x,y,z);}
        }
        vtkfile.end_vertices();

        vtkfile.begin_cells();
        for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != Tri.finite_cells_end(); ++cell) {
		bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal() && cell->vertex(2)->info().isReal()  && cell->vertex(3)->info().isReal();
        	if (isDrawable){vtkfile.write_cell(cell->vertex(0)->info().id()-firstReal, cell->vertex(1)->info().id()-firstReal, cell->vertex(2)->info().id()-firstReal, cell->vertex(3)->info().id()-firstReal);}
        }
        vtkfile.end_cells();

	if (permeabilityMap){
	vtkfile.begin_data("Permeability",CELL_DATA,SCALARS,FLOAT);
	for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != Tri.finite_cells_end(); ++cell) {
		bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal() && cell->vertex(2)->info().isReal()  && cell->vertex(3)->info().isReal();
		if (isDrawable){vtkfile.write_data(cell->info().s);}
	}
	vtkfile.end_data();}
	else{
	vtkfile.begin_data("Pressure",CELL_DATA,SCALARS,FLOAT);
	for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != Tri.finite_cells_end(); ++cell) {
		bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal() && cell->vertex(2)->info().isReal()  && cell->vertex(3)->info().isReal();
		if (isDrawable){vtkfile.write_data(cell->info().p());}
	}
	vtkfile.end_data();}

	if (1){
	averageRelativeCellVelocity();
	vtkfile.begin_data("Velocity",CELL_DATA,VECTORS,FLOAT);
	for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != Tri.finite_cells_end(); ++cell) {
		bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal() && cell->vertex(2)->info().isReal()  && cell->vertex(3)->info().isReal();
		if (isDrawable){vtkfile.write_data(cell->info().averageVelocity()[0],cell->info().averageVelocity()[1],cell->info().averageVelocity()[2]);}
	}
	vtkfile.end_data();}
// 	/// Check this one, cell info()->cracked not defined yet
// 	if(1){
// // 	trickPermeability();
// 	vtkfile.begin_data("fracturedCells",CELL_DATA,SCALARS,FLOAT);
// 	for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != Tri.finite_cells_end(); ++cell) {
// 		bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal() && cell->vertex(2)->info().isReal()  && cell->vertex(3)->info().isReal();
// 		if (isDrawable){vtkfile.write_data(cell->info().crack);}
// 	}
// 	vtkfile.end_data();}
  }
  // e.g. vtk recorders
};


// Definition of the DFNFLowEngine
typedef TemplateFlowEngine<DFNCellInfo,DFNVertexInfo,DFNTesselation, DFNBoundingSphere> DFNFlowEngineT;
REGISTER_SERIALIZABLE(DFNFlowEngineT);
YADE_PLUGIN((DFNFlowEngineT));

class DFNFlowEngine : public DFNFlowEngineT
{
	public :
	  
	void trickPermeability();
	void trickPermeability (RTriangulation::Facet_circulator& facet,Real crackPermeability);
	void trickPermeability (RTriangulation::Finite_edges_iterator& edge,Real crackPermeability);

// Compiles and works without the doublication of the following functions:
// 	void action();
// 	void buildTriangulation (double pZero, Solver& flow);
// 	void buildTriangulation (Solver& flow);
// Up to here

	/// Old
// 	void checkForCracks( Solver& flow );
	/// Up to here

	YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(DFNFlowEngine,DFNFlowEngineT,"This is an enhancement of the FlowEngine for intact and fractured rocks that takes into acount pre-existing discontinuities and bond breakage between particles and multiplies the local conductivity around the broken link",
	((Real, newAttribute, 0,,"useless example, Input values to be implemented: newCrackResidualPermeability, jointResidualPermeability, crackOpeningFactor"))
	,
	,
	,
	.def("trickPermeability",&DFNFlowEngineT::trickPermeability,"measure the mean trickPermeability in the period")
// 	.def("clearImposedPressure",&TemplateFlowEngine::clearImposedPressure,"Clear the list of points with pressure imposed.")
	)
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(DFNFlowEngine);
YADE_PLUGIN((DFNFlowEngine));


/// Redifinition of buildTriangulation()
// Emergency Change
// void DFNFlowEngine::action()
// {
//        if ( !isActivated ) return;
//         timingDeltas->start();
// 	setPositionsBuffer(true);
// 	timingDeltas->checkpoint ( "Position buffer" );
//         if (first) {
// 	  if (multithread) setPositionsBuffer(false);
// 	  buildTriangulation(pZero,*solver);
// 	  initializeVolumes(*solver);
// 	  backgroundSolver=solver;
// 	  backgroundCompleted=true;
// 	}
// 	solver->ompThreads = ompThreads>0? ompThreads : omp_get_max_threads();
// 
//         timingDeltas->checkpoint ( "Triangulating" );
// 	updateVolumes ( *solver );
//         timingDeltas->checkpoint ( "Update_Volumes" );
// 	
//         epsVolCumulative += epsVolMax;
// 	retriangulationLastIter++;
// 	if (!updateTriangulation) updateTriangulation = // If not already set true by another function of by the user, check conditions
// 		(defTolerance>0 && epsVolCumulative > defTolerance) || retriangulationLastIter>meshUpdateInterval;
// 
//         ///compute flow and and forces here
// 	if (pressureForce){
// 		solver->gaussSeidel(scene->dt);
// 		timingDeltas->checkpoint ( "Gauss-Seidel (includes matrix construct and factorization in single-thread mode)" );
// 		solver->computeFacetForcesWithCache();}
//         timingDeltas->checkpoint ( "compute_Forces" );
//         ///Application of vicscous forces
//         scene->forces.sync();
// 	timingDeltas->checkpoint ( "forces.sync()" );
// 	computeViscousForces ( *solver );
// 	timingDeltas->checkpoint ( "viscous forces" );
// 	Vector3r force;
// 	Vector3r torque;
//         FiniteVerticesIterator verticesEnd = solver->T[solver->currentTes].Triangulation().finite_vertices_end();
//         for ( FiniteVerticesIterator vIt = solver->T[solver->currentTes].Triangulation().finite_vertices_begin(); vIt !=  verticesEnd; vIt++ ) {
// 		force = pressureForce ? Vector3r ( vIt->info().forces[0],vIt->info().forces[1],vIt->info().forces[2] ): Vector3r(0,0,0);
// 		torque = Vector3r(0,0,0);
//                 if (shearLubrication || viscousShear){
// 			force = force + solver->shearLubricationForces[vIt->info().id()];
// 			torque = torque + solver->shearLubricationTorques[vIt->info().id()];
// 			if (pumpTorque)
// 				torque = torque + solver->pumpLubricationTorques[vIt->info().id()];
// 		}
// 		if (twistTorque)
// 			torque = torque + solver->twistLubricationTorques[vIt->info().id()];
// 		if (normalLubrication)
// 			force = force + solver-> normalLubricationForce[vIt->info().id()];
// 		scene->forces.addForce ( vIt->info().id(), force);
// 		scene->forces.addTorque ( vIt->info().id(), torque);
//         }
//         ///End compute flow and forces
//         timingDeltas->checkpoint ( "Applying Forces" );
// 	int sleeping = 0;
// 	if (multithread && !first) {
// 		while (updateTriangulation && !backgroundCompleted) { /*cout<<"sleeping..."<<sleeping++<<endl;*/
// 		  sleeping++;
// 		boost::this_thread::sleep(boost::posix_time::microseconds(1000));}
// 		if (debug && sleeping) cerr<<"sleeping..."<<sleeping<<endl;
// 		if (updateTriangulation || (ellapsedIter>(0.5*meshUpdateInterval) && backgroundCompleted)) {
// 			if (debug) cerr<<"switch flow solver"<<endl;
// 			if (useSolver==0) LOG_ERROR("background calculations not available for Gauss-Seidel");
// 			if (fluidBulkModulus>0) solver->interpolate (solver->T[solver->currentTes], backgroundSolver->T[backgroundSolver->currentTes]);
// 			solver=backgroundSolver;
// 			backgroundSolver = shared_ptr<FlowSolver> (new FlowSolver);
// 			if (metisForced) {backgroundSolver->eSolver.cholmod().nmethods=1; backgroundSolver->eSolver.cholmod().method[0].ordering=CHOLMOD_METIS;}
// 			//Copy imposed pressures/flow from the old solver
// 			backgroundSolver->imposedP = vector<pair<CGT::Point,Real> >(solver->imposedP);
// 			backgroundSolver->imposedF = vector<pair<CGT::Point,Real> >(solver->imposedF);
// 			if (debug) cerr<<"switched"<<endl;
// 			setPositionsBuffer(false);//set "parallel" buffer for background calculation 
// 			backgroundCompleted=false;
// 			retriangulationLastIter=ellapsedIter;
// 			updateTriangulation=false;
// 			epsVolCumulative=0;
// 			ellapsedIter=0;
// 			boost::thread workerThread(&TemplateFlowEngine::backgroundAction,this);
// 			workerThread.detach();
// 			if (debug) cerr<<"backgrounded"<<endl;
// 			initializeVolumes(*solver);
// 			computeViscousForces(*solver);
// 			if (debug) cerr<<"volumes initialized"<<endl;
// 		}
// 		else {
// 			if (debug && !backgroundCompleted) cerr<<"still computing solver in the background, ellapsedIter="<<ellapsedIter<<endl;
// 			ellapsedIter++;
// 		}
// 	} else {
// 	        if (updateTriangulation && !first) {
// 			buildTriangulation (pZero, *solver);
// 			initializeVolumes(*solver);
// 			computeViscousForces(*solver);
//                		updateTriangulation = false;
// 			epsVolCumulative=0;
// 			retriangulationLastIter=0;
// 			ReTrg++;}
//         }
//         first=false;
//         timingDeltas->checkpoint ( "triangulate + init volumes" );
// }
// 
// 
// void DFNFlowEngine::buildTriangulation ( Solver& flow )
// {
//         buildTriangulation ( 0.f,flow );
// }
// void DFNFlowEngine::buildTriangulation ( double pZero, Solver& flow )
// {
//  	if (first) flow.currentTes=0;
//         else {
//                 flow.currentTes=!flow.currentTes;
//                 if (debug) cout << "--------RETRIANGULATION-----------" << endl;
//         }
// 	flow.resetNetwork();
// 	initSolver(flow);
// 
//         addBoundary ( flow );
//         triangulate ( flow );
//         if ( debug ) cout << endl << "Tesselating------" << endl << endl;
//         flow.T[flow.currentTes].compute();
// 
//         flow.defineFictiousCells();
// 	// For faster loops on cells define this vector
// 	flow.T[flow.currentTes].cellHandles.clear();
// 	flow.T[flow.currentTes].cellHandles.reserve(flow.T[flow.currentTes].Triangulation().number_of_finite_cells());
// 	FiniteCellsIterator cell_end = flow.T[flow.currentTes].Triangulation().finite_cells_end();
// 	int k=0;
// 	for ( FiniteCellsIterator cell = flow.T[flow.currentTes].Triangulation().finite_cells_begin(); cell != cell_end; cell++ ){
// 		flow.T[flow.currentTes].cellHandles.push_back(cell);
// 		cell->info().id=k++;}//define unique numbering now, corresponds to position in cellHandles
//         flow.displayStatistics ();
//         flow.computePermeability();
// 	//This virtual function does nothing yet, derived class may overload it to make permeability different (see DFN engine)
// 	cout << " Starting trickPermeability () function from the copied function" << endl;
// 	this->trickPermeability();
// 	cout << " Just finished trickPermeability() function from the copied function" << endl;
//         porosity = flow.vPoralPorosity/flow.vTotalPorosity;
// 
//         boundaryConditions ( flow );
//         flow.initializePressure ( pZero );
// 	
//         if ( !first && !multithread && (useSolver==0 || fluidBulkModulus>0)) flow.interpolate ( flow.T[!flow.currentTes], flow.T[flow.currentTes] );
//         if ( waveAction ) flow.applySinusoidalPressure ( flow.T[flow.currentTes].Triangulation(), sineMagnitude, sineAverage, 30 );
//         if (normalLubrication || shearLubrication || viscousShear) flow.computeEdgesSurfaces();
// }

/// Old
// void DFNFlowEngine::trickPermeability()
// {
// // 	if ( debug ) cout << "Checking for cracked Edges/Interactions " << endl;
// // 	cout << "Entering checkForCracks() || Checking for cracked Edges/Interactions " << endl; // debug comment
// 	Solver& flow = *solver;
// 	int nbCrackedCells=0;
// 	const JCFpmPhys* jcfpmphys;
// 	const shared_ptr<InteractionContainer> interactions = scene->interactions;
// 	const RTriangulation& Tri = solver->T[solver->currentTes].Triangulation();
// 	Real newCrackPermMultiplier=1000;
// 	for ( FiniteCellsIterator cell = flow.T[flow.currentTes].Triangulation().finite_cells_begin(); cell != flow.T[flow.currentTes].Triangulation().finite_cells_end(); cell++ ) {
// 		if ( cell->info().isFictious ) continue;
// 		int bondedEdges = 0;
// 		for ( int g=0;g<3;g++ ) {
// 			for ( int h=g+1;h<4;h++ ) {
// 				const shared_ptr<Interaction>& interaction=interactions->find( cell->vertex(g)->info().id(), cell->vertex(h)->info().id() );
// 				if ( (interaction!=0) && (interaction->isReal()) ) {
// 					jcfpmphys = YADE_CAST<JCFpmPhys*>(interaction->phys.get());
// 					if (jcfpmphys->isCohesive) bondedEdges+=1;
// 				}
// 			}
// 		}
//         	if ( bondedEdges < 4 ) { cell->info().crack=1; nbCrackedCells+=1; }
// 	}
// 	for ( FiniteCellsIterator cell = flow.T[flow.currentTes].Triangulation().finite_cells_begin(); cell != flow.T[flow.currentTes].Triangulation().finite_cells_end(); cell++ ) {
// 		if ( cell->info().isFictious ) continue;
// 		for (int i=0 ; i<4; i++) {
// 		  if ( cell->info().crack>0 && cell->neighbor(i)->info().crack>0 ) { (cell->info().kNorm())[i]*=newCrackPermMultiplier; (cell->neighbor(i)->info().kNorm())[Tri.mirror_index(cell, i)]*=newCrackPermMultiplier;
// 		  }
// 		}
// 	}
// 	flow.computePermeability(); // needed to update local permeabilities
// 	cout << "Leaving checkForCracks() now.. " << endl;; // debug comment
// }
/// Up to here

// The trickPermeability function exist already in the FlowEngine but it's empty. Here whe can define what it will do.
  void DFNFlowEngine::trickPermeability (RTriangulation::Facet_circulator& facet, Real crackPermeability)
  {
  	const RTriangulation& Tri = solver->T[solver->currentTes].Triangulation();
  	const CellHandle& cell1 = facet->first;
  	const CellHandle& cell2 = facet->first->neighbor(facet->second);
  	if ( Tri.is_infinite(cell1) || Tri.is_infinite(cell2)) cerr<<"Infinite cell found in trickPermeability, should be handled somehow, maybe"<<endl;
// 	continue; ???
  	cell1->info().kNorm()[facet->second] = crackPermeability;
  	cell2->info().kNorm()[Tri.mirror_index(cell1, facet->second)] = crackPermeability;
	cout<< "Permeability changed to ...some facets ------------------------------" << endl;
  }

//Circulate over the facets on which, a specific edge is an instance.
void DFNFlowEngine::trickPermeability (RTriangulation::Finite_edges_iterator& edge, Real crackPermeability)
{
	const RTriangulation& Tri = solver->T[solver->currentTes].Triangulation();
	RTriangulation::Facet_circulator facet1 = Tri.incident_facets(*edge);
	RTriangulation::Facet_circulator facet0=facet1++;
	trickPermeability (facet0,crackPermeability);
	while ( facet1!=facet0 ) {trickPermeability(facet1,crackPermeability); facet1++;}
	cout << "circulator finished ------------------------------- " << endl;
}


void DFNFlowEngine::trickPermeability()
{
	Real crackPermeability=0.0000001; // This is the value for the permeability multiplier of the newly created cracks
// 	Real preExistingJointPermMultiplier = __ ; // This is the value for the permeability multiplier of the pre-existing cracks

	const RTriangulation& Tri = solver->T[solver->currentTes].Triangulation();

// 	#ifdef GS_OPEN_MP
// 	pragma openmp shared(edge,..) private() schedule(static??, chunck)
//     #endif
	
	/// Loop logic:
 	/// Loop into edges 
	/// Find the corresponding interaction using find() function
	///  If the interaction is cracked or injoint
	///  Apply trickPermeability() == Circulate around the broken link and change the conductivity of all facets on which this link is an instance
	
	cout << "Checking for cracked Edges/Interactions " << endl; // Debuging comment Delete afterwards 
	
	const JCFpmPhys* jcfpmphys;
	const shared_ptr<InteractionContainer> interactions = scene->interactions;
	FiniteEdgesIterator edge = Tri.finite_edges_begin();
	int TestOverInteractionsLooped= 0; // Debuging comment Delete afterwards 
	for( ; edge!= Tri.finite_edges_end(); ++edge) {
		const VertexInfo& vi1=(edge->first)->vertex(edge->second)->info();
		const VertexInfo& vi2=(edge->first)->vertex(edge->third)->info();
		const shared_ptr<Interaction>& interaction=interactions->find( vi1.id(),vi2.id() );
		if (interaction && interaction->isReal()) {
			jcfpmphys = YADE_CAST<JCFpmPhys*>(interaction->phys.get());
			if ( (!jcfpmphys->isCohesive) ) {trickPermeability(edge,crackPermeability);}; // changed reco
			cout << " permeability enhancement around edge between spheres: || " << vi1.id() << " and " << vi2.id() << endl;// Debuging comment Delete afterwards  
		}
		TestOverInteractionsLooped +=1; // Debuging comment Delete afterwards 
		cout << " number of interactions Looped up to now || " << TestOverInteractionsLooped << endl; // Debuging comment Delete afterwards  
	}
}

#endif //DFNFLOW
#endif //FLOWENGINE

/// Cell-based approach
// FlowBoundingSphere.ipp, computePermeability() ... 
// if ( cell->info().crack>0 && neighbour_cell->info().crack>0 ) { (cell->info().k_norm())[j]*=kCrackFactor; (neighbour_cell->info().k_norm())[Tri.mirror_index(cell, j)]*=kCrackFactor; }

// FlowEngine.cpp, checkForCracks() ... 


/// Some instructions..
	// vertex 1 of edge e is:
	// e->first->vertex(e->second)
	// vertex 2:
	// e->first->vertex(e->third)
	// better:
	// (ed_it->first)->vertex(ed_it->second)->info() /// from computeEdgesSurfaces in flowBoundingSphere.ipp
	
