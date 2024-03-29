/******************************************************************************/
/*!
\file   contraints.graph.h
\author Dmitri Volper
\par    email: dvolper\@digipen.edu
\date   Wed 15 Apr 2009 12:19:18 PM PDT
\brief  
  Class representing data structure containing constraints and variables
  for Constraint Satisfaction Problem.

  Implements several types of methods
  1) insertion of constraints and variables
  2) preprocessing constraints and variables for fast retrieval
  3) methods to retrieve sub-collections of constraints and/or variables.
  4) checks
*/
/******************************************************************************/

/******************************************************************************/
/*!
  \class Variable
  \brief  
  Class representing all constraints for Constraint Satisfaction Problem.
  Implements:
	insertion of variables
  	insertion of constraints
	retrieval of all constraints by pointer to a Variable
	retrieval of all Variable connected to a given Variable
*/
/******************************************************************************/
#if defined(_MSC_VER)
// warning C4290: C++ exception specification ignored except to indicate
// a function is not __declspec(nothrow)
# pragma warning(disable:4290)
//warning C4996: 'std::xxx' was declared deprecated # pragma warning(disable:4996)
#endif 


#ifndef CONSTRAINT_GRAPH_H
#define CONSTRAINT_GRAPH_H
#include <vector>
#include <map>
#include <set>
#include <deque>
#include <iostream>
#include <string>
#include <cstdarg>
#include <cassert>
#include <utility>
#include <algorithm>
#include "contraints.h"


//constraint graph - used in CSP Problem
//just a data structure to simplify access to constraints from variables
//and vice versa
template <typename T>
class ConstraintGraph {
	public:
		//typedef's
		typedef T Constraint;
		typedef typename Constraint::Variable Variable;

		//ctors,dtors
		ConstraintGraph();
		~ConstraintGraph();

		//insertion methods
		////////////////////////////////////////////////////////////
		void InsertVariable( Variable& var );
		////////////////////////////////////////////////////////////
		void InsertConstraint( const Constraint & c ) throw (const char *);

		////////////////////////////////////////////////////////////
		//pre-build collections
		void PreProcess();
		
		//retrieval methods
		////////////////////////////////////////////////////////////
		//set of Variables which are connected to a given Variable 
		//by a constraint 
		const typename std::set<Variable*> & GetNeighbors( Variable* p_var );
		////////////////////////////////////////////////////////////
		//all constraints that use the given variables
		const typename std::vector<const Constraint*> & GetConstraints( Variable* p_var ) const throw (const char *);
		////////////////////////////////////////////////////////////
		//vector of Constraints connecting 2 given Variables
		const typename std::set<const Constraint*> & GetConnectingConstraints( Variable* p_var1, Variable* p_var2 );
		////////////////////////////////////////////////////////////
		//vector of all Variables
		const typename std::vector<Variable*> & GetAllVariables( ) const;

		//checks
		////////////////////////////////////////////////////////////
		//activate/disactivate constraints
		void CheckActivity();
		////////////////////////////////////////////////////////////
		//detect solution (that is -- all variable are assigned)
		bool AllVariablesAssigned() const;

		////////////////////////////////////////////////////////////
		void Print() const;
		////////////////////////////////////////////////////////////
	private:
		//data-structures to simplify life

		//mapping from variable to constraints the variable is used in.
		//that is constraints that depend on the variable
		// (x, {1,2,3,4,5})
		std::map<Variable*,std::vector<const Constraint*> > var2constr;
		//all variables
		std::vector<Variable*> vars;
		//all constraints
		std::vector<Constraint*> constraints;//mainly for print

		//mapping from variable x to a set of variables Z
		//each z in Z is connected to x by a constraint
		//that is: there is a constraint that uses both variables x and z 
		std::map<Variable*, std::set<Variable*> > neighbors;

		//for each pair of variables all constraints that use both of them
		//used in consistency checking
		std::map< std::pair<Variable*,Variable*>, std::set<const Constraint*> > connecting_constraints;
		//for internal use only
		std::map<std::string,Variable*> name2vars;
};

#ifdef INLINE_CONSTRAINT_GRAPH
//#warning "INFO - inlining ConstraintGraph methods"
#define INLINE inline
#else   
//#warning "INFO - NOT inlining ConstraintGraph methods"
#define INLINE 
#endif

template <typename T>
ConstraintGraph<T>::ConstraintGraph() :
	var2constr(),
	vars(),
	constraints(),
	neighbors(),
	connecting_constraints(),
	name2vars()
{
}

template <typename T>
ConstraintGraph<T>::~ConstraintGraph() {
	typename std::vector<Constraint*>::const_iterator
		b_constr = constraints.begin();
	typename std::vector<Constraint*>::const_iterator
		e_constr = constraints.end();

	for (;b_constr != e_constr;++b_constr) {
		delete *b_constr;
	}

}
////////////////////////////////////////////////////////////
//activate/disactivate constraints
template <typename T>
INLINE
void ConstraintGraph<T>::CheckActivity() {
	typename std::vector<Constraint*>::iterator
		b_constr = constraints.begin();
	typename std::vector<Constraint*>::iterator
		e_constr = constraints.end();

	for (;b_constr != e_constr;++b_constr) {
		(*b_constr)->SetActive();
	}
}
////////////////////////////////////////////////////////////
//build the above data-structures
template <typename T>
void ConstraintGraph<T>::PreProcess() {
	typename std::vector<Variable*>::const_iterator b_vars = vars.begin();
	typename std::vector<Variable*>::const_iterator e_vars = vars.end();
	for (;b_vars != e_vars;++b_vars) {
		const std::vector<const Constraint*> & constr = GetConstraints(*b_vars);

		//for all constraints current variable is used in
		//find all other variables and put them in the set
		//set ensures no duplicates
		std::set<Variable*> neigh;
		typename std::vector<const Constraint*>::const_iterator b_constr = constr.begin();
		typename std::vector<const Constraint*>::const_iterator e_constr = constr.end();

		for (;b_constr != e_constr;++b_constr) {
			typename std::vector< Variable*>::const_iterator b_vars2 = (*b_constr)->GetVars().begin();
			typename std::vector< Variable*>::const_iterator e_vars2 = (*b_constr)->GetVars().end();
			for (;b_vars2 != e_vars2;++b_vars2) {
				if (*b_vars2 != *b_vars)//excluding the given variable
				{
					neigh.insert(*b_vars2);
					//current constraint "connects" 2 variables
					std::pair<Variable*, Variable*> varPair(*b_vars, *b_vars2);
					connecting_constraints[varPair].insert(*b_constr);
				}

			}
		}
		neighbors[*b_vars] = neigh;
	}
}

////////////////////////////////////////////////////////////
//set of Variables which are connected to a given Variable 
//by a constraint 
template <typename T>
INLINE
const typename std::set<typename ConstraintGraph<T>::Variable*>&
ConstraintGraph<T>::GetNeighbors(typename ConstraintGraph<T>::Variable* p_var) {
	return neighbors[p_var];
}
////////////////////////////////////////////////////////////
//set of unassigned Variables which are connected to a given Variable 
//by a constraint 
//template <typename T>
//unsigned ConstraintGraph<T>::GetCountUnassignedNeighbors( Variable* p_var ) {
//	unsigned result;
//	typename std::set<Variable*>::const_iterator b_n = neighbors[p_var].begin();
//	typename std::set<Variable*>::const_iterator e_n = neighbors[p_var].end();
//	for ( ; b_n != e_n; ++b_n ) {
//		result += ( (*b_n)->IsAssigned() ? 0:1 );
//	}
//	return result;
//}
////////////////////////////////////////////////////////////
//vector of Constraints connecting 2 given Variables
template <typename T>
INLINE
const typename std::set<const typename ConstraintGraph<T>::Constraint*>&
ConstraintGraph<T>::GetConnectingConstraints(
	typename ConstraintGraph<T>::Variable* p_var1,
	typename ConstraintGraph<T>::Variable* p_var2)
{
	std::pair<Variable*, Variable*> varPair(p_var1, p_var2);
	return connecting_constraints[varPair];
}
////////////////////////////////////////////////////////////
//detect dead-end
//template <typename T>
//bool ConstraintGraph<T>::CheckPossibleAssignments() const {
//	typename std::vector<typename ConstraintGraph<T>::Variable*>::const_iterator 
//		b_vars = vars.begin();
//	typename std::vector<typename ConstraintGraph<T>::Variable*>::const_iterator 
//		e_vars = vars.end();
//	for ( ;b_vars!=e_vars;++b_vars) { 
//		if ( ! (*b_vars)->IsImpossible() ) return false;
//	}
//	return true;
//}
////////////////////////////////////////////////////////////
template <typename T>
void ConstraintGraph<T>::InsertVariable(typename ConstraintGraph<T>::Variable& var) {
	//Variable* p_var = new Variable(var);
	Variable* p_var;
	p_var = &var;

	vars.push_back(p_var);
	name2vars[var.Name()] = p_var;
	std::pair<Variable*, std::vector<const Constraint*> > pair(p_var, std::vector<const Constraint*>());
	var2constr.insert(pair); //initially empty 

							 //std::cout << __FILE__ << " " << __LINE__ << " var address " << p_var << std::endl;
							 //std::cout << *p_var << std::endl; 
							 //std::cout << std::endl;
}
////////////////////////////////////////////////////////////
template <typename T>
void ConstraintGraph<T>::InsertConstraint(const Constraint & c) throw (const char *) {
	Constraint* p_c = c.clone();
	//std::cout << "local constraint " << *p_c << std::endl;
	const std::vector<Variable*> & vars_in_constraint = p_c->GetVars();
	//check we know all variables
	typename std::vector<Variable*>::const_iterator b = vars_in_constraint.begin();
	typename std::vector<Variable*>::const_iterator e = vars_in_constraint.end();

	//insert constraint as an outgoing to all variables
	//used in the constraint
	for (;b != e;++b) {
		typename std::map<std::string, Variable*>::iterator it = name2vars.find((*b)->Name());
		if (it == name2vars.end()) {
			std::cout << "unknown variable name " << (*b)->Name() << std::endl;
			throw "constraint uses unknown variable name";
		}
		//successfully found variable
		var2constr[it->second].push_back(p_c);
	}
	constraints.push_back(p_c);
}
////////////////////////////////////////////////////////////
template <typename T>
INLINE
const typename std::vector<const typename ConstraintGraph<T>::Constraint*>&
ConstraintGraph<T>::GetConstraints(typename ConstraintGraph<T>::Variable* p_var) const
throw (const char *)
{
	typename
		std::map< Variable*, std::vector<const Constraint*> >::const_iterator it = var2constr.find(p_var);
	if (it != var2constr.end()) return it->second;
	else throw "cannot find variable in the graph";
}
////////////////////////////////////////////////////////////
//vector of all Variables
template <typename T>
INLINE
const typename std::vector<typename ConstraintGraph<T>::Variable*>&
ConstraintGraph<T>::GetAllVariables() const {
	return vars;
}
////////////////////////////////////////////////////////////
//detect solution (that is -- all variables are assigned)
template <typename T>
INLINE
bool ConstraintGraph<T>::AllVariablesAssigned() const {
	typename std::vector<Variable*>::const_iterator
		b_vars = vars.begin();
	typename std::vector<Variable*>::const_iterator
		e_vars = vars.end();
	for (;b_vars != e_vars;++b_vars) {
		if (!(*b_vars)->IsAssigned()) return false;
	}
	return true;
}
////////////////////////////////////////////////////////////
template <typename T>
void ConstraintGraph<T>::Print() const {
	typename std::vector<Variable*>::const_iterator
		b_vars = vars.begin();
	typename std::vector<Variable*>::const_iterator
		e_vars = vars.end();
	std::cout << "+------------------+\n";
	std::cout << "|  ConstraintGraph |\n";
	std::cout << "+------------------+\n";
	std::cout << "|-->Variables :\n";
	for (;b_vars != e_vars;++b_vars) {
		std::cout << "|   " << **b_vars << std::endl;
		//std::cout << "   "; (*b_vars)->Print();
	}

	typename std::vector<Constraint*>::const_iterator
		b_constr = constraints.begin();
	typename std::vector<Constraint*>::const_iterator
		e_constr = constraints.end();
	std::cout << "|-->Constraints :\n";
	for (;b_constr != e_constr;++b_constr) {
		std::cout << "|   "; std::cout << **b_constr << std::endl;
	}
	std::cout << "+------------------+\n";
}
////////////////////////////////////////////////////////////
#undef INLINE


#endif
