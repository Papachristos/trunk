 
/*************************************************************************
*  Copyright (C) 2014 by Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>     *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/
	//TODO: Cubic Law on fractures
	//TODO: Add a " CrackWidth" factor (I guess, of the same order as the permeabFactor (??) )
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
// 	trickPermeability();
// 	vtkfile.begin_data("fracturedCells",CELL_DATA,SCALARS,FLOAT);
// 	for (Finite_cells_iterator cell = Tri.finite_cells_begin(); cell != Tri.finite_cells_end(); ++cell) {
// 		bool isDrawable = cell->info().isReal() && cell->vertex(0)->info().isReal() && cell->vertex(1)->info().isReal() && cell->vertex(2)->info().isReal()  && cell->vertex(3)->info().isReal();
// 		if (isDrawable){vtkfile.write_data(cell->info().crack);}
// 	}
// 	vtkfile.end.data();}
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
	void trickPermeability (RTriangulation::Facet_circulator& facet,Real newCrackPermMultiplier);
	void trickPermeability (RTriangulation::Finite_edges_iterator& edge,Real newCrackPermMultiplier);

	YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(DFNFlowEngine,DFNFlowEngineT,"This is an enhancement of the FlowEngine that takes into acount pre-existing discontinuities and bond breakage between particles and multiplies the local conductivity around the broken link",
	((Real, myNewAttribute, 0,,"useless example, Input values to be implemented: newCrackPermMultiplier, jointCrackPermMultiplie, residualCrackPermMultiplier, residualJointPermMultiplier, crackOpeningFactor"))
	,
	,
	,
	//.def("trickPermeability",&DFNFlowEngine::trickPermeability,"measure the mean trickPermeability in the period")
	)
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(DFNFlowEngine);
YADE_PLUGIN((DFNFlowEngine));

// The trickPermeability function exist already in the FlowEngine but it's empty. Here whe can define what it will do.
void DFNFlowEngine::trickPermeability (RTriangulation::Facet_circulator& facet, Real newCrackPermMultiplier)
{
	const RTriangulation& Tri = solver->T[solver->currentTes].Triangulation();
	const CellHandle& cell1 = facet->first;
	const CellHandle& cell2 = facet->first->neighbor(facet->second);
	if ( Tri.is_infinite(cell1) || Tri.is_infinite(cell2)) cerr<<"Infinite cell found in trickPermeability, should be handled somehow, maybe"<<endl;
	cell1->info().kNorm()[facet->second] = newCrackPermMultiplier;
	cell2->info().kNorm()[Tri.mirror_index(cell1, facet->second)] = newCrackPermMultiplier;
}

// Circulate over the facets on which, a specific edge is an instance.
void DFNFlowEngine::trickPermeability(RTriangulation::Finite_edges_iterator& edge, Real newCrackPermMultiplier)
{
	const RTriangulation& Tri = solver->T[solver->currentTes].Triangulation();
	RTriangulation::Facet_circulator facet1 = Tri.incident_facets(*edge);
	RTriangulation::Facet_circulator facet0=facet1++;
	trickPermeability(facet0,newCrackPermMultiplier);
	while ( facet1!=facet0 ) {trickPermeability(facet1,newCrackPermMultiplier); facet1++;}
}


void DFNFlowEngine::trickPermeability()
{

	Real newCrackPermMultiplier=10000; // This is the value for the permeability multiplier of the newly created cracks
// 	Real preExistingJointPermMultiplier = __ ; // This is the value for the permeability multiplier of the pre-existing cracks

	const RTriangulation& Tri = solver->T[solver->currentTes].Triangulation();

// 	#ifdef GS_OPEN_MP
// 	pragma openmp shared(edge,..) private() schedule(static??, chunck)
//     #endif
	
	/// Loop logic:
 	// Loop into edges 
	// Find the corresponding interaction using find() function
	//  If the interaction is cracked or injoint
	//  Apply trickPermeability() == Circulate around the broken link and change the conductivity of all facets on which this link is an instance
	
	cout << "Checking for cracked Edges/Interactions " << endl;
	
	const JCFpmPhys* jcfpmphys;
	const shared_ptr<InteractionContainer> interactions = scene->interactions;
	FiniteEdgesIterator edge = Tri.finite_edges_begin();
	for( ; edge!= Tri.finite_edges_end(); ++edge) {
		const VertexInfo& vi1=(edge->first)->vertex(edge->second)->info();
		const VertexInfo& vi2=(edge->first)->vertex(edge->third)->info();
		const shared_ptr<Interaction>& interaction=interactions->find( vi1.id(),vi2.id() );
		if (interaction->isReal() ) {
			jcfpmphys = YADE_CAST<JCFpmPhys*>(interaction->phys.get());
			if ( (!jcfpmphys->isCohesive) ) {trickPermeability(edge,newCrackPermMultiplier);}; // changed reco
			cout << " permeability enhancement around edge: || " << endl;
		}
	}
}

#endif //DFNFLOW
#endif //FLOWENGINE

// Some instructions..
	// vertex 1 of edge e is:
	// e->first->vertex(e->second)
	// vertex 2:
	// e->first->vertex(e->third)
	// better:
	// (ed_it->first)->vertex(ed_it->second)->info() /// from computeEdgesSurfaces in flowBoundingSphere.ipp
	
