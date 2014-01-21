// © 2004 Janek Kozicki <cosurgi@berlios.de>
// © 2007 Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>
// © 2008 Václav Šmilauer <eudoxos@arcig.cz>

#include"Ig2_Sphere_Sphere_ScGeom.hpp"
#include<yade/pkg/dem/ScGeom.hpp>
#include<yade/pkg/common/Sphere.hpp>
#include<yade/core/Scene.hpp>
#include<yade/lib/base/Math.hpp>
#include<yade/core/Omega.hpp>
#include<yade/pkg/common/InteractionLoop.hpp>

bool Ig2_Sphere_Sphere_ScGeom::go(	const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c)
{
	TIMING_DELTAS_START();
	const Se3r& se31=state1.se3; const Se3r& se32=state2.se3;
	const Sphere *s1=static_cast<Sphere*>(cm1.get()), *s2=static_cast<Sphere*>(cm2.get());
	Vector3r normal=(se32.position+shift2)-se31.position;
	if (!c->isReal() && !force) {//don't fast-check distance if geometry will be updated anyway
		Real penetrationDepthSq=pow(interactionDetectionFactor*(s1->radius+s2->radius),2) - normal.squaredNorm();
		if (penetrationDepthSq<0) {
			TIMING_DELTAS_CHECKPOINT("Ig2_Sphere_Sphere_ScGeom");
			return false;
		}
	}
	shared_ptr<ScGeom> scm;
	bool isNew = !c->geom;
	if(!isNew) scm=YADE_PTR_CAST<ScGeom>(c->geom);
	else { scm=shared_ptr<ScGeom>(new ScGeom()); c->geom=scm; scm->penetrationInc=1.;}
	Real norm=normal.norm(); normal/=norm; // normal is unit vector now
#ifdef YADE_DEBUG
	if(norm==0) throw runtime_error(("Zero distance between spheres #"+lexical_cast<string>(c->getId1())+" and #"+lexical_cast<string>(c->getId2())+".").c_str());
#endif
	Real penetrationDepth=s1->radius+s2->radius-norm;
	if (interactionDetectionFactor!=1){
		scm->contactPoint=se31.position+(s1->radius-(s1->radius/(s1->radius+s2->radius))*penetrationDepth)*normal;//OC-p*un;
	}
	else scm->contactPoint=se31.position+(s1->radius-0.5*penetrationDepth)*normal;//0.5*(pt1+pt2);
	#ifdef NORATCH2
	scm->alpha = (s1->radius+s2->radius)/norm;
// 	scm->penetrationInc=penetrationDepth-scm->penetrationDepth;
// 	scm->objective = objective;
	#endif
	scm->penetrationDepth=penetrationDepth;
	scm->radius1=s1->radius;
	scm->radius2=s2->radius;
	scm->shearScheme=shearScheme;
// 	scm->qBased = qBased;
	scm->precompute(state1,state2,scene,c,normal,isNew,shift2,avoidGranularRatcheting);
	TIMING_DELTAS_CHECKPOINT("Ig2_Sphere_Sphere_ScGeom");
	return true;
}

bool Ig2_Sphere_Sphere_ScGeom::goReverse(	const shared_ptr<Shape>& cm1,
								const shared_ptr<Shape>& cm2,
								const State& state1,
								const State& state2,
								const Vector3r& shift2,
								const bool& force,
								const shared_ptr<Interaction>& c)
{
	return go(cm1,cm2,state2,state1,-shift2,force,c);
}

YADE_PLUGIN((Ig2_Sphere_Sphere_ScGeom));

bool Ig2_Sphere_Sphere_ScGeom6D::go( const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c)
{
	bool isNew = !c->geom;
	if (Ig2_Sphere_Sphere_ScGeom::go(cm1,cm2,state1,state2,shift2,force,c)){//the 3 DOFS from ScGeom are updated here
 		if (isNew) {//generate a 6DOF interaction from the 3DOF one generated by Ig2_Sphere_Sphere_ScGeom
			shared_ptr<ScGeom6D> sc (new ScGeom6D());
			*(YADE_PTR_CAST<ScGeom>(sc)) = *(YADE_PTR_CAST<ScGeom>(c->geom));
			c->geom=sc; sc->penetrationInc=1.;}
		if (updateRotations) YADE_PTR_CAST<ScGeom6D>(c->geom)->precomputeRotations(state1,state2,isNew,creep);
		return true;
	}
	else return false;
}

bool Ig2_Sphere_Sphere_ScGeom6D::goReverse( const shared_ptr<Shape>& cm1, const shared_ptr<Shape>& cm2, const State& state1, const State& state2, const Vector3r& shift2, const bool& force, const shared_ptr<Interaction>& c)
{
	return go(cm1,cm2,state2,state1,-shift2,force,c);
}

YADE_PLUGIN((Ig2_Sphere_Sphere_ScGeom6D));
