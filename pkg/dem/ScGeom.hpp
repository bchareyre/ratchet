// © 2004 Olivier Galizzi <olivier.galizzi@imag.fr>
// © 2008 Václav Šmilauer <eudoxos@arcig.cz>
// © 2006 Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>

#pragma once
#include<yade/core/Interaction.hpp>
#include<yade/core/IGeom.hpp>
#include<yade/core/State.hpp>
#include<yade/lib/base/Math.hpp>
#include<yade/pkg/dem/DemXDofGeom.hpp>

/*! Class representing geometry of two bodies in contact.
 *
 */

#define NORATCH2

class ScGeom: public GenericSpheresContact {
	private:
		//cached values
		Vector3r twist_axis;//rotation vector around normal
		Vector3r orthonormal_axis;//rotation vector in contact plane
		enum { OBJECTIVE=1, NORATCHETING=2, MID_STEP_NORMAL=4, EXACT_SHEAR=8, EXACT_ROTATE=16 };
	public:
		// inherited from GenericSpheresContact: Vector3r& normal;
		Real &radius1, &radius2;
		virtual ~ScGeom();
		inline ScGeom& operator= (const ScGeom& source){
			normal=source.normal; contactPoint=source.contactPoint;
			twist_axis=source.twist_axis; orthonormal_axis=source.orthonormal_axis;
			radius1=source.radius1; radius2=source.radius2;
			penetrationDepth=source.penetrationDepth; shearInc=source.shearInc;
			return *this;}

		//!precompute values of shear increment and interaction rotation data. Update contact normal to the currentNormal value. Precondition : the value of normal is not updated outside (and before) this function.
		void precompute(const State& rbp1, const State& rbp2, const Scene* scene, const shared_ptr<Interaction>& c, const Vector3r& currentNormal, bool isNew, const Vector3r& shift2, bool avoidGranularRatcheting=true);

		//! Rotates a "shear" vector to keep track of contact orientation. Returns reference of the updated vector.
		Vector3r& rotate(Vector3r& tangentVector) const;
		const Vector3r& shearIncrement() const {return shearInc;}

		// Add method which returns the relative velocity (then, inside the contact law, this can be split into shear and normal component). Handle periodicity.
		Vector3r getIncidentVel(const State* rbp1, const State* rbp2, Real dt, const Vector3r& shiftVel, const Vector3r& shift2, bool avoidGranularRatcheting=true);
		// Implement another version of getIncidentVel which does not handle periodicity.
		Vector3r getIncidentVel(const State* rbp1, const State* rbp2, Real dt, bool avoidGranularRatcheting=true);
		// Add function to get the relative angular velocity (useful to determine bending moment at the contact level)
		Vector3r getRelAngVel(const State* rbp1, const State* rbp2, Real dt);

		// return a vector orthogonal to argument
		Vector3r orthVect (const Vector3r& v);
		// rotate a vector using Rodriguez relation. In all cases, axisAngle is NOT a unit vector
		// passing more arguments favors performance, the functions are in the order of increasing speed
		Vector3r rrot (const Vector3r& axisAngle, const Vector3r& v) const;
		Vector3r rrot (const Vector3r& axisAngle, const Real& angle, const Vector3r& v) const;
		Vector3r rrot (const Vector3r& axisAngle, const Real& angle, const Real& sinA, const Real& cosA, const Vector3r& v) const;
		Vector3r rrotWithSin (const Vector3r& axisSinAngle, const Vector3r& v) const;
		Vector3r rrotWithSin (const Vector3r& axisSinAngle, const Real& sinA, const Vector3r& v) const;
		Vector3r rrotWithSin (const Vector3r& axisSinAngle, const Real& sinA, const Real& cosA, const Vector3r& v) const;
		
		
		// convenience version to be called from python
		Vector3r getIncidentVel_py(shared_ptr<Interaction> i, bool avoidGranularRatcheting);
		Vector3r getRelAngVel_py(shared_ptr<Interaction> i);

		YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(ScGeom,GenericSpheresContact,"Class representing :yref:`geometry<IGeom>` of a contact point between two :yref:`bodies<Body>`. It is more general than sphere-sphere contact even though it is primarily focused on spheres interactions (reason for the 'Sc' namming); it is also used for representing contacts of a :yref:`Sphere` with non-spherical bodies (:yref:`Facet`, :yref:`Plane`,  :yref:`Box`, :yref:`ChainedCylinder`), or between two non-spherical bodies (:yref:`ChainedCylinder`). The contact has 3 DOFs (normal and 2×shear) and uses incremental algorithm for updating shear.\n\nWe use symbols $\\vec{x}$, $\\vec{v}$, $\\vec{\\omega}$ respectively for position, linear and angular velocities (all in global coordinates) and $r$ for particles radii; subscripted with 1 or 2 to distinguish 2 spheres in contact. Then we define branch length and unit contact normal\n\n.. math::\n\n\tl=||\\vec{x}_2-\\vec{x}_1||, \\vec{n}=\\frac{\\vec{x}_2-\\vec{x}_1}{||\\vec{x}_2-\\vec{x}_1||}\n\nThe relative velocity of the spheres is then\n\n.. math::\n\n\t\\vec{v}_{12}=\\frac{r_1+r_2}{l}(\\vec{v}_2-\\vec{v}_1) -(r_2 \\vec{\\omega}_2 + r_1\\vec{\\omega}_1)\\times\\vec{n}\n\nwhere the fraction multplying translational velocities is to make the definition objective and avoid ratcheting effects (see :yref:`Ig2_Sphere_Sphere_ScGeom.avoidGranularRatcheting`). The shear component is\n\n.. math::\n\n\t\\vec{v}_{12}^s=\\vec{v}_{12}-(\\vec{n}\\cdot\\vec{v}_{12})\\vec{n}.\n\nTangential displacement increment over last step then reads\n\n.. math::\n\n\t\\Delta\\vec{x}_{12}^s=\\Delta t \\vec{v}_{12}^s.",
		((Real,penetrationDepth,0,(Attr::readonly),"Penetration distance of spheres (positive if overlapping)"))
		
		((int,shearScheme,3,,"Selectively define the algorithm to be used for the shear components.\n\n"
		"== ===============================================================\n"
		"1  Objective? uses the historical Cundall's expression with constant radius if 0 (discards options 2,4,8), '1' let other options be used.\n"
		"2  avoid granular ratcheting? (see :yref:`Ig2_Sphere_Sphere_ScGeom::avoidGranularRatcheting`) \n"
		"4  use mid-step normal instead of current step for the projection of relative velocity \n"
		"8  compute exact relative shear instead of projecting the shearVelocity*dt. \n"
		"16  use exact equations for updating contact orientation, gives more precision but takes more time.…\n"
		"== ===============================================================\n\n"
		"By default, the mask is 3, using a fast non-ratcheting algorithm, at the price of precision. Other meaningful values are:\n"
		"0, for experimenting the historical formulation \n"
		"1, for the most commonly used equations \n"
		"5, like above but with mid-step normal (as in PFC) \n"
		"7, non-ratcheting with mid-step normal \n"
		"11, exact non-ratcheting formulation, approximated contact rotation \n"
		"16,17,19,21,23,27, exact contact rotation combined with the above. Exact rotations are not really useful without exact shear displacements, since the result would be inexact anyway.\n"
		"Forbidden/inconsistent: 2,4,8,9,12,14."
		))
		#ifdef NORATCH2
		((Real,penetrationInc,0,(Attr::noSave|Attr::readonly),"Increment of normal displacement between two steps. Used by non-ractcheting code."))
		((Vector3r,midStepN,Vector3r::Zero(),," - "))
		((Vector3r,us,Vector3r::Zero(),," - "))
		((Vector3r,t1,Vector3r::Zero(),," - "))
		((Real,alpha,1,," - "))
		#endif
		((Vector3r,shearInc,Vector3r::Zero(),(Attr::noSave|Attr::readonly),"Shear displacement increment in the last step"))
		,
		/* extra initializers */ ((radius1,GenericSpheresContact::refR1)) ((radius2,GenericSpheresContact::refR2)),
		/* ctor */ createIndex(); twist_axis=orthonormal_axis=Vector3r::Zero();,
		/* py */ .def("incidentVel",&ScGeom::getIncidentVel_py,(python::arg("i"),python::arg("avoidGranularRatcheting")=true),"Return incident velocity of the interaction (see also :yref:`Ig2_Sphere_Sphere_ScGeom.avoidGranularRatcheting` for explanation of the ratcheting argument).")
		.def("relAngVel",&ScGeom::getRelAngVel_py,(python::arg("i")),"Return relative angular velocity of the interaction.")
	);
	REGISTER_CLASS_INDEX(ScGeom,GenericSpheresContact);
};
REGISTER_SERIALIZABLE(ScGeom);

class ScGeom6D: public ScGeom {
	public:
		virtual ~ScGeom6D();
		const Real& getTwist() const {return twist;}
		const Vector3r& getBending() const {return bending;}
		void precomputeRotations(const State& rbp1, const State& rbp2, bool isNew, bool creep=false);
		void initRotations(const State& rbp1, const State& rbp2);
		YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(ScGeom6D,ScGeom,"Class representing :yref:`geometry<IGeom>` of two :yref:`bodies<Body>` in contact. The contact has 6 DOFs (normal, 2×shear, twist, 2xbending) and uses :yref:`ScGeom` incremental algorithm for updating shear.",
		((Quaternionr,initialOrientation1,Quaternionr(1.0,0.0,0.0,0.0),(Attr::readonly),"Orientation of body 1 one at initialisation time |yupdate|"))
		((Quaternionr,initialOrientation2,Quaternionr(1.0,0.0,0.0,0.0),(Attr::readonly),"Orientation of body 2 one at initialisation time |yupdate|"))
		((Quaternionr,twistCreep,Quaternionr(1.0,0.0,0.0,0.0),(Attr::readonly),"Stored creep, substracted from total relative rotation for computation of elastic moment |yupdate|"))
		((Real,twist,0,(Attr::noSave | Attr::readonly),"Elastic twist angle of the contact."))
		((Vector3r,bending,Vector3r::Zero(),(Attr::noSave | Attr::readonly),"Bending at contact as a vector defining axis of rotation and angle (angle=norm)."))
		,
		/* extra initializers */,
		/* ctor */ createIndex();,
 		/* py */
	);
	REGISTER_CLASS_INDEX(ScGeom6D,ScGeom);
};
REGISTER_SERIALIZABLE(ScGeom6D);


class ChCylGeom6D: public ScGeom6D {
	public:
		virtual ~ChCylGeom6D();
		State fictiousState1;
		State fictiousState2;
		Real relPos1, relPos2;
		YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(ChCylGeom6D,ScGeom6D,"Test",
		/*attributes*/
		,
		/* extra initializers */,
		/* ctor */ createIndex();,
		/* py */
	);
	REGISTER_CLASS_INDEX(ChCylGeom6D,ScGeom6D);
};
REGISTER_SERIALIZABLE(ChCylGeom6D);























