 
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
	//stupid change
};


typedef CGT::_Tesselation<CGT::TriangulationTypes<DFNVertexInfo,DFNCellInfo> > DFNTesselation;
// We add all the new/complementary things for FlowBoundingSphere.ipp in the DFNboundingSphere
// always declare public what we call from another class
class DFNBoundingSphere : public CGT::FlowBoundingSphere<DFNTesselation>
{
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

	YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(DFNFlowEngine,DFNFlowEngineT,"documentation here",
	((Real, myNewAttribute, 0,,"useless example"))
	,
// 	DFNFlowEngineT()
	,
	,
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

// Circulate over the facets on which, a specific edge is an incident.
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

	Real newCrackPermMultiplier=10; // This is the value for the permeability multiplier of the newly created cracks
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
	for(int edge; edge!= Tri.finite_edges_end(); ++edge) {
		const shared_ptr<Interaction>& interaction=interactions->find( edge->first->vertex(edge->second), edge->first->vertex(edge->third) );
		if (interaction->isReal() ) {
			jcfpmphys = YADE_CAST<JCFpmPhys*>(interaction->phys.get());
			if (jcfpmphys->isCohesive) || interaction == 0 ) {trickPermeability(edge,newCrackPermMultiplier)};
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
