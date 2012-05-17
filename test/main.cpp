#ifndef _MAIN_H_
#define _MAIN_H_

#include <iostream>
#include <fstream>
#include <stdlib.h>

#include <xperience/core>

xp::Error initTest(xp::Context& ctx);
xp::Error finishTest(xp::Context& ctx);
xp::Error testExperiment(xp::Context& ctx);

int main(int argc, char *argv[])
{
  xp::Scope scope("Test Agent", std::string("NO_SHA1"));

  scope.registerExperiment("testExperiment", &testExperiment, &initTest, &finishTest);


  xp::Error err=scope.parseCommandLine(argc, argv);

  return err;
}

xp::Error initTest(xp::Context& ctx)
{
  return xp::kNoErr;
}

xp::Error finishTest(xp::Context& ctx)
{
  return xp::kNoErr;
}

xp::Error testExperiment(xp::Context& ctx)
{
  CTX_LOG(ctx, info) << mulog::format::red << "Started Experiment";

  mulog::logger l = ctx.openDataLogger("log");

  LLOG(l,info) << mulog::format::green << "Test";

  return xp::kNoErr;
}

#endif /* _MAIN_H_ */
