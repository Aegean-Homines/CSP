#ifndef CSP_H
#define CSP_H
#include <vector>
#include <set>
#include <iostream>
#include <string>
#include <cstdarg>
#include <cassert>
#include <deque>
#include <utility>
#include <algorithm>
#include <limits>
#include <map>
#include <ctime>

template <typename C>
struct Arc {
	typedef typename C::Variable Variable;
	//non-owning semantics, therefore DO NOT need big-4
	Variable* x;
	Variable* y;
	const C* c;
	Arc(Variable* x,Variable* y,const C* c) : x(x),y(y),c(c) {}
	bool operator< (const Arc<C>& rhs) const {
		if ( x<rhs.x) return true;
		if ( x==rhs.x && y<rhs.y) return true;
		if ( x==rhs.x && y==rhs.y && c<rhs.c) return true;

		return false;
	}
};

template <typename T> 
class CSP {
		//typedef's for intenal use
		typedef typename T::Constraint      Constraint;
		typedef typename T::Variable        Variable;
		typedef typename T::Variable::Value Value;
	public:
		////////////////////////////////////////////////////////////
		//counters
		////////////////////////////////////////////////////////////
		//ctor
		CSP(T &cg);

		//get the number of found solutions
		int GetSolutionCounter() const { return solution_counter; }
		//get the number of recursive calls - for debugging
		int GetRecursiveCallCounter() const { return recursive_call_counter; }
		//get the number of variable assigns in Solve* - for debugging
		int GetIterationCounter() const { return iteration_counter; }

		//CSP counting
		bool SolveFC_count(unsigned level);
		//CSP solver, brute force - no forward checking
		bool SolveDFS(unsigned level);
		//CSP solver, uses forward checking
		bool SolveFC(unsigned level);
		//CSP solver, uses arc consistency
		bool SolveARC(unsigned level);
	private:
		//2 versions of forward checking algorithms
		bool ForwardChecking(Variable *x);
		//load states (available values) of all unassigned variables 
		void LoadState(std::map<Variable*, std::set<typename CSP<T>::Variable::Value> >& saved) 
			//note: "CSP<T>::" in "std::set<typename CSP<T>::Variable::Value>" required for MSC
		const ;
		//save states (available values) of all unassigned variables 
		//except the current
		std::map<Variable*, std::set<typename Variable::Value> > 
			SaveState(Variable* x) const;
		//check the current (incomplete) assignment for satisfiability
		bool AssignmentIsConsistent( Variable* p_var ) const;
		//insert pair 
		//(neighbors of the current variable, the current variable)
		//for all y~x insert (y,x)
		//into arc-consistency queue
		void InsertAllArcsTo( Variable* cv );

		//AIMA p.146 AC-3 algorithm
		bool CheckArcConsistency(Variable* x);
		//CHECK that for each value of x there is a value of y 
		//which makes all constraints involving x and y satisfiable
		bool RemoveInconsistentValues(Variable* x,Variable* y,const Constraint* c);
		//choose next variable for assignment
		//choose the one with minimum remaining values
		Variable* MinRemValByDomain();
		Variable* MinRemValByConstraint();
		Variable* FirstAvailable();

		//choose next variable for assignment
		//choose the one with max degree
		Variable* MaxDegreeHeuristic();

		// EGEMEN
		// For debugging
		void Print() const {
			const std::vector<Variable*>& all_vars = cg.GetAllVariables();
			typename std::vector<Variable*>::const_iterator b_vars = all_vars.begin();
			typename std::vector<Variable*>::const_iterator e_vars = all_vars.end();
			std::cout << "------------------------" << std::endl;
			std::cout << "@@@@ SOLUTION FOUND @@@@" << std::endl;
			for (; b_vars != e_vars; ++b_vars) {
				(*b_vars)->Print();
				std::cout << std::endl;
			}
		}

		//data
		//deque of arcs (2 Variables connected through a Constraint)
		std::set< Arc<Constraint> > arc_consistency;
		T &cg;
		int solution_counter,recursive_call_counter,iteration_counter;
};


#ifdef INLINE_CSP
//#warning "INFO - inlining CSP methods"
#define INLINE inline
#else   
//#warning "INFO - NOT inlining CSP methods"
#define INLINE 
#endif

////////////////////////////////////////////////////////////
//CSP constructor
template <typename T>
CSP<T>::CSP(T &cg) :
	arc_consistency(),
	cg(cg),
	solution_counter(0),
	recursive_call_counter(0),
	iteration_counter(0)
{
}

////////////////////////////////////////////////////////////
//CSP solver, brute force - no forward checking
template <typename T>
bool CSP<T>::SolveDFS(unsigned level) {
	++recursive_call_counter;

	// This only happens when we call the solve for a non-existing level
	// That is this recursion is called by the last variable assigned
	// However, if all variables are assigned without a problem at this point
	// Then we have the solution of the problem with the assigned values
	if (cg.AllVariablesAssigned()) {
		return true;
	}
	

#ifdef DEBUG
	std::cout << "entering SolveDFS (level " << level << ")\n";
#endif
	Variable* var_to_assign = NULL;
	//choose a variable by MRV
	//Variable* var_to_assign = MaxDegreeHeuristic();
#ifndef FIRST_AVAILABLE

	var_to_assign = MinRemValByDomain();
#else
	var_to_assign = FirstAvailable();

	if (!var_to_assign) {
		return false;
	}
#endif

	// Get all the constraints that we need to satisfy
	const std::vector<const Constraint*> & allConstraints = cg.GetConstraints(var_to_assign);
	typename std::vector<const Constraint*>::const_iterator b_const = allConstraints.begin();
	typename std::vector<const Constraint*>::const_iterator e_const = allConstraints.end();
	// Get variables domain values (Save state)
	const std::set<Value>& var_domain = var_to_assign->GetDomain();
	typename std::set<Value> original_domain = var_domain;


	// While the possible values of domain set is not empty
	while (!var_to_assign->IsImpossible()) {
		++iteration_counter;

		// Assign the first element in the set to this variable
		var_to_assign->Assign();

		bool allSatisfiable = true;

		// Reset the iterator for the next loop
		b_const = allConstraints.begin();
		// Does this assignment satisfy all constraints
		for (; (b_const != e_const) && (allSatisfiable); ++b_const) {
			allSatisfiable &= (*b_const)->Satisfiable();
		}

		if (!allSatisfiable) {
			var_to_assign->RemoveValue(var_to_assign->GetValue()); //deletes from the original 
			var_to_assign->UnAssign();
			continue;
		}
		else if (!SolveDFS(level + 1)) {
			var_to_assign->RemoveValue(var_to_assign->GetValue()); //deletes from the original 
			var_to_assign->UnAssign();
			continue;
		}
		else {
#ifdef DEBUG
			Print();
#endif
			return true;
		}

		//This should've worked as well, not different from whats going on at the above code snippet
		/*if (!allSatisfiable || !SolveDFS(level + 1)) { //If this value doesn't break satisfiable condition and sub branches also work
			var_to_assign->RemoveValue(var_to_assign->GetValue()); //deletes from the original 
			var_to_assign->UnAssign();
		}
		else {
			Print(allVariables);
			return true;
		}*/

	}

	if (level == 0) { // The first variable took all the values and its children checked all possible values
		std::cout << "NO SOLUTION EXISTS!" << std::endl;
	}
	else { //Not first -> rollback domain state
		var_to_assign->SetDomain(original_domain);
	}
	return false;

}


////////////////////////////////////////////////////////////
//CSP solver, uses forward checking
template <typename T>
bool CSP<T>::SolveFC(unsigned level) {
	++recursive_call_counter;
#ifdef DEBUG
	std::cout << "entering SolveFC (level " << level << ")\n";
#endif

	if (cg.AllVariablesAssigned()) {
		return true;
	}

	Variable* var_to_assign = NULL;

#ifndef FIRST_AVAILABLE
	var_to_assign = MinRemValByDomain();
#else
	var_to_assign = FirstAvailable();
	if (!var_to_assign) {
		return false;
	}
#endif
	const std::vector<const Constraint*> & allConstraints = cg.GetConstraints(var_to_assign);
	typename std::vector<const Constraint*>::const_iterator b_const = allConstraints.begin();
	typename std::vector<const Constraint*>::const_iterator e_const = allConstraints.end();
	// Get variables domain values (Save state)
	const std::set<Value>& var_domain = var_to_assign->GetDomain();
	typename std::set<Value> original_domain = var_domain;

	// While the possible values of domain set is not empty
	while (!var_to_assign->IsImpossible()) {
		++iteration_counter;

		// Assign the first element in the set to this variable
		var_to_assign->Assign();

		bool allSatisfiable = true;

		// Reset the iterator for the next loop
		b_const = allConstraints.begin();
		// Does this assignment satisfy all constraints
		for (; (b_const != e_const) && (allSatisfiable); ++b_const) {
			allSatisfiable &= (*b_const)->Satisfiable();
		}

		if (!allSatisfiable) {
			var_to_assign->RemoveValue(var_to_assign->GetValue()); //deletes from the original 
			var_to_assign->UnAssign();
			continue;
		}

		std::map< typename CSP<T>::Variable*, std::set<typename CSP<T>::Variable::Value> > savedState;
		savedState = SaveState(var_to_assign); //save the current state

		bool isImpossible = ForwardChecking(var_to_assign);

		if (!isImpossible && SolveFC(level + 1)) {
			return true;
		}
		else {
			LoadState(savedState);
			var_to_assign->RemoveValue(var_to_assign->GetValue()); //deletes from the original 
			var_to_assign->UnAssign();
			continue;
		}

		/*if (!isImpossible) { // No solution with this variable
			LoadState(savedState);
			var_to_assign->RemoveValue(var_to_assign->GetValue()); //deletes from the original 
			var_to_assign->UnAssign();
			continue;
		}
		else if(!SolveFC(level + 1)) {
			LoadState(savedState);
			var_to_assign->RemoveValue(var_to_assign->GetValue()); //deletes from the original 
			var_to_assign->UnAssign();
			continue;
		}
		else {
#ifdef DEBUG
			Print();
#endif
			return true;
		}*/
		
	}

	if (level == 0) { // The first variable took all the values and its children checked all possible values
		std::cout << "NO SOLUTION EXISTS!" << std::endl;
	}

	var_to_assign->SetDomain(original_domain);

	return false;

}
////////////////////////////////////////////////////////////
//CSP solver, uses arc consistency
template <typename T>
bool CSP<T>::SolveARC(unsigned level) {
	++recursive_call_counter;
	std::cout << "entering SolveARC (level " << level << ")\n";




	/*
	//choose a variable by MRV
	Variable* var_to_assign = MinRemVal();





	loop( ... ) {
	++iteration_counter;



	}
	*/

	return false;


}


template <typename T>
INLINE
bool CSP<T>::ForwardChecking(Variable *x) {
	// Get the neighbors of current variable
	const std::set<Variable*> & neighbors = cg.GetNeighbors(x);
	typename std::set<Variable*>::const_iterator b_neig = neighbors.begin();
	typename std::set<Variable*>::const_iterator e_neig = neighbors.end();

	while (b_neig != e_neig) {
		Variable * currentNeighbor = *(b_neig++);
		if (currentNeighbor->IsAssigned()) //pass assigned values
			continue;

		// Get all constraints between this and current neighbor
		const std::set<const Constraint*> & connecting_constrains = cg.GetConnectingConstraints(x, currentNeighbor);
		typename std::set<const Constraint*>::const_iterator b_const = connecting_constrains.begin();
		typename std::set<const Constraint*>::const_iterator e_const = connecting_constrains.end();

		// Get all the values of the neighbor
		const std::set<Value>& all_vals = currentNeighbor->GetDomain();
		typename std::set<Value>::const_iterator b_vals = all_vals.begin();
		typename std::set<Value>::const_iterator e_vals = all_vals.end();

		// Go through all the values
		while (b_vals != e_vals) {
			Value currentValue = *(b_vals++);
			// Assign each value to the neighbor
			currentNeighbor->Assign(currentValue);
			b_const = connecting_constrains.begin();

			// Loop through all the constraints for this value
			while (b_const != e_const) {
				Constraint const * currentConstraint = *(b_const++);

				// Does this value satisfy the constraint
				if (!currentConstraint->Satisfiable()) {
					// If not remove the value and set neighbor to its unassigned state
					currentNeighbor->RemoveValue(currentNeighbor->GetValue());
					currentNeighbor->UnAssign(); //This shouldn't be necessary but I don't wanna deal with possible bugs
					if (currentNeighbor->SizeDomain() == 0) { //If that was the last value then everything is wrong
						return true;
					}
					break;
				}

			}
		}

		if (currentNeighbor->IsAssigned())
			currentNeighbor->UnAssign();

	}


	return false;
}
////////////////////////////////////////////////////////////
//load states (available values) of all unassigned variables 
template <typename T>
void CSP<T>::LoadState(std::map<Variable*, std::set<typename CSP<T>::Variable::Value> >& saved) const
{
	typename std::map<Variable*, std::set<typename Variable::Value> >::iterator b_result = saved.begin();
	typename std::map<Variable*, std::set<typename Variable::Value> >::iterator e_result = saved.end();

	for (; b_result != e_result; ++b_result) {
		//std::cout << "loading state for " 
		//<< b_result->first->Name() << std::endl;
		(*b_result).first->SetDomain((*b_result).second);
	}
}


////////////////////////////////////////////////////////////
//save states (available values) of all unassigned variables 
//except the current
template <typename T>
INLINE
std::map< typename CSP<T>::Variable*, std::set<typename CSP<T>::Variable::Value> >
CSP<T>::SaveState(typename CSP<T>::Variable* x) const {
	std::map<Variable*, std::set<typename Variable::Value> > result;

	const std::vector<Variable*>& all_vars = cg.GetAllVariables();
	typename std::vector<Variable*>::const_iterator b_all_vars = all_vars.begin();
	typename std::vector<Variable*>::const_iterator e_all_vars = all_vars.end();
	for (; b_all_vars != e_all_vars; ++b_all_vars) {
		if (!(*b_all_vars)->IsAssigned() && *b_all_vars != x) {
			//std::cout << "saving state for " 
			//<< (*b_all_vars)->Name() << std::endl;
			result[*b_all_vars] = (*b_all_vars)->GetDomain();
		}
	}
	return result;
}
////////////////////////////////////////////////////////////
//check the current (incomplete) assignment for satisfiability
template <typename T>
INLINE
bool CSP<T>::AssignmentIsConsistent(Variable* p_var) const {

	return false;








}
////////////////////////////////////////////////////////////
//insert pair 
//(neighbors of the current variable, the current variable)
//current variable is th variable that just lost some values
// for all y~x insert (y,x)
//into arc-consistency queue
template <typename T>
INLINE
void CSP<T>::InsertAllArcsTo(Variable* cv) {











}
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
//AIMA p.146 AC-3 algorithm
template <typename T>
INLINE
bool CSP<T>::CheckArcConsistency(Variable* x) {


	return false;









}
////////////////////////////////////////////////////////////
//CHECK that for each value of x there is a value of y 
//which makes all constraints involving x and y satisfiable
template <typename T>
INLINE
bool CSP<T>::RemoveInconsistentValues(Variable* x, Variable* y, const Constraint* c) {

	return false;










}
////////////////////////////////////////////////////////////
//choose next variable for assignment
//choose the one with minimum remaining values
template <typename T>
INLINE
typename CSP<T>::Variable* CSP<T>::MinRemValByDomain() {
	// Get all the variables and loop through them
	const std::vector<Variable*> & allVariables = cg.GetAllVariables();
	typename std::vector<Variable*>::const_iterator b_vars = allVariables.begin();
	typename std::vector<Variable*>::const_iterator e_vars = allVariables.end();
	Variable* var_min = NULL;
	unsigned int max_size = std::numeric_limits<unsigned int>::max();

	for (; b_vars != e_vars; ++b_vars) {
		Variable* var_temp = *b_vars;
		const std::set<Value>& var_domain = var_temp->GetDomain();
		// If the variable in question is not assigned and if it's domain is smaller than the last max we found pick that one
		if (!var_temp->IsAssigned() && var_domain.size() < max_size) {
			var_min = var_temp;
			max_size = var_domain.size();
		}
	}
	return var_min;

}

template <typename T>
INLINE
typename CSP<T>::Variable* CSP<T>::MinRemValByConstraint() {
	// Get all the variables and loop through them
	const std::vector<Variable*> & allVariables = cg.GetAllVariables();
	typename std::vector<Variable*>::const_iterator b_vars = allVariables.begin();
	typename std::vector<Variable*>::const_iterator e_vars = allVariables.end();
	Variable* var_min = NULL;
	unsigned int min_size = std::numeric_limits<unsigned int>::min();

	for (; b_vars != e_vars; ++b_vars) {
		Variable* var_temp = *b_vars;
		const std::vector<const Constraint*> & allConstraints = cg.GetConstraints(var_temp);
		// If the variable in question is not assigned and if it's domain is smaller than the last max we found pick that one
		if (!var_temp->IsAssigned() && allConstraints.size() > min_size) {
			var_min = var_temp;
		}
	}
	return var_min;

}

template<typename T>
INLINE
typename CSP<T>::Variable* CSP<T>::FirstAvailable()
{
	const std::vector<Variable*> & allVariables = cg.GetAllVariables();
	typename std::vector<Variable*>::const_iterator b_vars = allVariables.begin();
	typename std::vector<Variable*>::const_iterator e_vars = allVariables.end();
	Variable* var_to_assign = NULL;
	// Get the first unassigned variable for this level
	for (; b_vars != e_vars; ++b_vars) {
		if (!(*b_vars)->IsAssigned()) {
			var_to_assign = *b_vars;
			break;
		}
	}

	return var_to_assign;
}


////////////////////////////////////////////////////////////
//choose next variable for assignment
//choose the one with max degree
template <typename T>
typename CSP<T>::Variable* CSP<T>::MaxDegreeHeuristic() {

	return NULL;












}
#undef INLINE


#endif
