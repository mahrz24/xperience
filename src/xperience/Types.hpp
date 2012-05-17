#ifndef _TYPES_H_
#define _TYPES_H_

#include <boost/function.hpp>

namespace xp
{
  class Context;

typedef enum ParamType
{
  kUnknown = 0,
  kDouble,
  kInt,
  kString,
  kBool
} ParamType;

typedef enum Error
{
kNoErr = 0,
  kNotInitialized = -1,
  kAlreadyInitialized = 1,
  kNameAlreadyRegistered,
  kExperimentNotFound,
  kNameNotFound,
  kFileNotFound,
  kSetNotFound,
  kParamFileFormatError,
  kParamParseError,
  kSweepRangeNotSpecified,
  kFileSystemError,
  kExperimentError,
  kExperimentInitError,
  kExperimentFinishError

  } Error;

typedef enum Phase
{
  kInit = 0,
  kRunning,
  kFinish
} Phase;

typedef boost::function<Error (Context&)> Callback;

}

#endif /* _TYPES_H_ */
