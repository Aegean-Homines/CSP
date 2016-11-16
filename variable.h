/******************************************************************************/
/*!
\file   variable.h
\author Dmitri Volper
\par    email: dvolper\@digipen.edu
\date   Tue 31 Mar 2009 01:51:20 PM PDT
\brief  
  Class representing a variable for Constraint Satisfaction Problem.
  Implements domain, assigned/not assigned state.
  Variable has a name and a unique id (id is not currently used)
  
*/
/******************************************************************************/

/******************************************************************************/
/*!
  \class Variable
  \brief  
  Class representing a variable for Constraint Satisfaction Problem.
  Implements domain, assigned/not assigned state.
  Variable has a name and a unique id (id is not currently used)

    Operations include:

*/
/******************************************************************************/
#if defined(_MSC_VER)
// warning C4290: C++ exception specification ignored except to indicate
// a function is not __declspec(nothrow)
# pragma warning(disable:4290)
//warning C4996: 'std::xxx' was declared deprecated # pragma warning(disable:4996)
#endif 


#ifndef VARIABLE_H
#define VARIABLE_H
#include <vector>
#include <set>
#include <fstream>
#include <string>
#include <iostream>
#include <algorithm>

class VariableException : public std::exception {
	std::string msg;
	public:
	VariableException(const std::string & _msg) : msg(_msg) {}
	const char * what() const throw () { return msg.c_str(); }
	virtual ~VariableException() throw () {}
};

class Variable {
	public:
		//! variable has to specify it's value type - used by other classes 
		typedef int Value;

	private:
		static unsigned nextid;

		//! symbol/name of the variable 
		std::string name;

		//! all remaining available values
		std::set<Value> domain;

		//! currently assigned value (invalid if is_assigned == false)
		Value assigned_value; 

		//! validity of the assigned_value 
		bool is_assigned; 

		//! unique ID (not used)
		unsigned id;


	public:
		Variable(const std::string & name, const std::vector<Value> & av):
			name(name),
			domain(av.begin(), av.end()),
			assigned_value(Value()),
			is_assigned(false),
			id(++nextid) {

		}
		friend std::ostream& operator<<(std::ostream& os, const Variable& v) {
			std::set<Variable::Value>::const_iterator b = v.GetDomain().begin();
			std::set<Variable::Value>::const_iterator e = v.GetDomain().end();
			os << "Variable \"" << v.Name() << "\" available values: ";
			for (;b != e;++b) { os << *b << " "; }
			os << std::endl;
			if (v.IsAssigned())
				os << " \nassigned value " << v.GetValue() << std::endl;
			return os;
		}
		const std::string & Name() const;
		void  RemoveValue(Value val) throw (VariableException);
		void  SetDomain(const std::set<Value>& vals);
		int   SizeDomain() const;
		const std::set<Value>& GetDomain() const;
		bool  IsImpossible() const;
		unsigned ID() const;
		bool  IsAssigned() const;
		void  Assign(Value val) throw (VariableException);
		void  Assign() throw (VariableException);
		void  UnAssign() throw (VariableException);
		Value GetMinValue() const throw (VariableException);
		Value GetMaxValue() const throw (VariableException);
		Value GetValue() const throw (VariableException);
		void  Print() const {
			std::cout << *this;
		}
};

#ifndef _WIN32
unsigned Variable::nextid = 0;
#endif

#ifdef INLINE_VARIABLE
	//#warning "INFO - inlining Variable methods"
	#include "variable.inl"
#endif

#ifndef INLINE_VARIABLE
	//#warning "INFO - NOT inlining Variable"
#include "variable.inl"
#endif

#endif
