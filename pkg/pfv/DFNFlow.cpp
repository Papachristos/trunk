 
/*************************************************************************
*  Copyright (C) 2014 by Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>     *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

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
	if(1){
// // 	trickPermeability();
	vtkfile.begin_data("fracturedCells",CELL_DATA,SCALARS,FLOAT);
	for (FiniteCellsIterator cell = Tri.finite_cells_begin(); cell != Tri.finite_cells_end(); ++cell) {
		bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal() && cell->vertex(2)->info().isReal()  && cell->vertex(3)->info().isReal();
		if (isDrawable){vtkfile.write_data(cell->info().crack);}
	}
	vtkfile.end_data();}
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
	void trickPermeability (RTriangulation::Facet_circulator& facet,Real crackPermeability,Real aperture, Real residualAperture);
	void trickPermeability (RTriangulation::Finite_edges_iterator& edge,Real crackPermeability,Real aperture, Real residualAperture);
	void setPositionsBuffer(bool current);

	YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(DFNFlowEngine,DFNFlowEngineT,"This is an enhancement of the FlowEngine for intact and fractured rocks that takes into acount pre-existing discontinuities and bond breakage between particles and multiplies the local conductivity around the broken link",
	((Real, newAttribute, 0,,"useless example, Input values to be implemented: newCrackResidualPermeability, jointResidualPermeability, crackOpeningFactor"))
	((bool, updatePositions,false,,"update particles positions when rebuilding the mesh (experimental)"))
	((Real, jointResidual, 0,,"Calibration parameter for residual aperture of joints"))
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
//In this version, we never update positions when !updatePositions, i.e. keep triangulating the same positions
void DFNFlowEngine::setPositionsBuffer(bool current)
{
	vector<posData>& buffer = current? positionBufferCurrent : positionBufferParallel;
	if (!updatePositions && buffer.size()>0) return;
	buffer.clear();
	buffer.resize(scene->bodies->size());
	shared_ptr<Sphere> sph ( new Sphere );
        const int Sph_Index = sph->getClassIndexStatic();
	FOREACH ( const shared_ptr<Body>& b, *scene->bodies ) {
                if (!b || ignoredBody==b->getId()) continue;
                posData& dat = buffer[b->getId()];
		dat.id=b->getId();
		dat.pos=b->state->pos;
		dat.isSphere= (b->shape->getClassIndex() ==  Sph_Index);
		if (dat.isSphere) dat.radius = YADE_CAST<Sphere*>(b->shape.get())->radius;
		dat.exists=true;
	}
}

// The trickPermeability function exist already in the FlowEngine but it's empty. Here whe can define what it will do.
  void DFNFlowEngine::trickPermeability (RTriangulation::Facet_circulator& facet, Real crackPermeability, Real aperture, Real residualAperture)
  {
  	const RTriangulation& Tri = solver->T[solver->currentTes].Triangulation();
  	const CellHandle& cell1 = facet->first;
  	const CellHandle& cell2 = facet->first->neighbor(facet->second);
  	if ( Tri.is_infinite(cell1) || Tri.is_infinite(cell2)) cerr<<"Infinite cell found in trickPermeability, should be handled somehow, maybe"<<endl;
// 	continue; ???
// 	cell1->info().kNorm()[facet->second] = crackPermeability;
  	cell1->info().kNorm()[facet->second] = pow((aperture+residualAperture),3) / (12 * viscosity); /// Timos - Check If we need to add density and viscosity - 
	if (cell1->info().kNorm()[Tri.mirror_index(cell1, facet->second)] !=8.33333e-08 ) {
	    cout << "Permeability set to : -- " << cell1->info().kNorm()[facet->second] << endl;
	}
//   	cell2->info().kNorm()[Tri.mirror_index(cell1, facet->second)] = crackPermeability;
	cell2->info().kNorm()[Tri.mirror_index(cell1, facet->second)] =pow((aperture+residualAperture),3) / (12 * viscosity);
	/// Check 
	cell1->info().crack= 1;
	cell2->info().crack= 1;
// 	cout<< "Permeability changed to ...some facets ------------------------------" << endl;
///  Cubic law:
// 
// 	dDem = jcfpmphys->dilation
//   	d = do a mean value
// 	dmax = for open fractures
// 	dmin = for closed fractures
// 	useCalibratedAperture = boolean
// 	apertureThreshold = a maximum value for aperture of open fractures
// 	
// if facet cracked {
// 
///   calculation of d
// 	if (useCalibratedAperture) {
// 	d = dDem
// 	}
//   	else if (d_dem<=0){
// 	 	d=dmin // "aperture" for closed fractures
// 	}
// 	else if (d_dem >= apertureThreshold) {
//   		dmax // higher allowed aperture value for open fractures
// 	}
//	else {
//   		d= do // a mean value for open fractures wiith aperture lower than apertureThreshold
// 	}
// 
///   calculation of beta
// 
// 
///   Caluclation of k
//	kNorm = ((beta * d).pow(3) )/(12 * viscosity)
// 
///   Mass flux
//	q = (d.pow(3)*(density*gravity))/12*viscosity *(Wij(Pi-Pj)/Lij) // Wij the width of the plates
// }
// 	
  }

//Circulate over the facets on which, a specific edge is an instance.
void DFNFlowEngine::trickPermeability (RTriangulation::Finite_edges_iterator& edge, Real crackPermeability, Real aperture, Real residualAperture )
{
	const RTriangulation& Tri = solver->T[solver->currentTes].Triangulation();
	RTriangulation::Facet_circulator facet1 = Tri.incident_facets(*edge);
	RTriangulation::Facet_circulator facet0=facet1++;
	trickPermeability (facet0,crackPermeability, aperture,residualAperture);
	while ( facet1!=facet0 ) {trickPermeability(facet1,crackPermeability, aperture, residualAperture); facet1++;}
// 	cout << "circulator finished ------------------------------- " << endl;
}


void DFNFlowEngine::trickPermeability()
{
	Real crackPermeability=0.0000001; // This is the value for the permeability multiplier of the newly created cracks
// 	Real preExistingJointPermMultiplier = __ ; // This is the value for the permeability multiplier of the pre-existing cracks

	const RTriangulation& Tri = solver->T[solver->currentTes].Triangulation();

// 	#ifdef GS_OPEN_MP
// 	pragma openmp shared(edge,..) private() schedule(static??, chunck)
//     #endif
	
// 	cout << "Checking for cracked Edges/Interactions " << endl; // Debuging comment Delete afterwards 
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
			if ( (jcfpmphys->interactionIsCracked || jcfpmphys->isOnJoint) ) {Real aperture=jcfpmphys->dilation; Real residualAperture = jcfpmphys->isOnJoint? jointResidual : 0; trickPermeability(edge,crackPermeability,aperture, residualAperture);}; // changed reco
// 			cout << " permeability enhancement around edge between spheres: || " << vi1.id() << " and " << vi2.id() << endl;// Debuging comment Delete afterwards  
		}
		TestOverInteractionsLooped +=1; // Debuging comment Delete afterwards 
// 		cout << " number of interactions Looped up to now || " << TestOverInteractionsLooped << endl; // Debuging comment Delete afterwards  
	}
}

#endif //DFNFLOW
#endif //FLOWENGINE

/// -------------------------------------------------------------------------------------------------------------------------------------------------------


/// Cell-based approach
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



/// Some instructions..
	// vertex 1 of edge e is:
	// e->first->vertex(e->second)
	// vertex 2:
	// e->first->vertex(e->third)
	// better:
	// (ed_it->first)->vertex(ed_it->second)->info() /// from computeEdgesSurfaces in flowBoundingSphere.ipp
	
