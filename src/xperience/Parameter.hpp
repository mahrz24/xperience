#ifndef _PARAMETER_H_
#define _PARAMETER_H_

#define PREFIX(p,s) (std::string(p) + s)

#include <tclap/CmdLine.h>
#include "Types.hpp"

namespace xp
{
  template<typename Type>
  struct TypeDescr
  {
    static const ParamType type = kUnknown;
    static const std::string descr() { return "argument"; };
    static const Type init() { return Type(); };
  };

  template<>
  struct TypeDescr<double>
  {
    static const ParamType type = kDouble;
    static const std::string descr() { return "float"; };
    static const double init() { return 0; };
  };

  template<>
  struct TypeDescr<int>
  {
    static const ParamType type = kInt;
    static const std::string descr() { return "integer"; };
    static const int init() { return 0; };
  };

  template<>
  struct TypeDescr<bool>
  {
    static const ParamType type = kBool;
    static const std::string descr() { return "bool"; };
    static const bool init() { return false; };
  };

  template<>
  struct TypeDescr<std::string>
  {
    static const ParamType type = kString;
    static const std::string descr() { return std::string("string"); };
    static const std::string init() { return std::string(""); };
  };

  class Parameter
  {
  public:
    Parameter() : type(kUnknown) {}
    Parameter(ParamType _type) : type(_type) {}

    ParamType type;

    virtual bool isSet() = 0;
    virtual void updateInternalValue() = 0;
  };

  template <typename Type>
  class TypedParameter : public Parameter
  {
  public:
    TypedParameter();

    TypedParameter(const std::string & parameterName,
              const std::string & parameterDescription);

    bool isSet() { return arg.isSet(); }
    void updateInternalValue() { val = arg.getValue(); }

    TCLAP::ValueArg<Type> arg;
    Type val;
  };

  template<typename Type>
  TypedParameter<Type>::TypedParameter() :
  Parameter(TypeDescr<Type>::type),
    arg("","" ,"", false,"", TypeDescr<Type>::descr()),
  val(TypeDescr<Type>::init())
  {}

  template<typename Type>
  TypedParameter<Type>::TypedParameter(const std::string & parameterName,
                                 const std::string & parameterDescription) :
                                 arg("", PREFIX("P_", parameterName) ,
                                     PREFIX("Experiment parameter: ", parameterDescription),
                                     false, TypeDescr<Type>::init() , TypeDescr<Type>::descr()),
      Parameter(TypeDescr<Type>::type),
      val(TypeDescr<Type>::init())
  {}
}

#endif /* _PARAMETER_H_ */
