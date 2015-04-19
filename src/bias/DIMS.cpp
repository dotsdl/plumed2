/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Bias.h"
#include "tools/Random.h"
#include "ActionRegister.h"
#include <ctime>
#include <cmath>
#include <random>

namespace PLMD{
namespace bias{

class DIMS : public Bias{

  // target value for CV
  std::vector<double> to;
  
  // current value for CV for this step
  std::vector<double> curr;

  // strength of the rejection potential
  std::vector<double> kappa;
 
  // width of the soft-ratcheting exponential distribution
  std::vector<double> phi;

  // seed for random number generator
  std::vector<int> seed;
public:
  DIMS(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(DIMS,"DIMS")

void DIMS::registerKeywords(Keywords& keys){
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","TO","array of target cutoffs for each CV");
  keys.add("compulsory","KAPPA","array of force constants for the rejection potential for each CV");
  keys.add("compulsory","PHI","array of soft-ratcheting parameters for each CV");
  keys.add("optional","SEED","array of seeds for the soft-ratcheting random variables");
  componentsAreNotOptional(keys);
  keys.addOutputComponent("accepted","default","whether the step was rejected or accepted");
  keys.addOutputComponent("_curr","default","one or multiple instances of this quantity will be refereceable elsewhere in the input file. "
                                 " These quantities will be named with the arguments of the bias followed by "
                                 "the character string _min. These quantities tell the user the minimum value assumed by rho_m(t).");
}

DIMS::DIMS(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao),
to(getNumberOfArguments(),0.0),
curr(getNumberOfArguments(),-1.0),
kappa(getNumberOfArguments(),0.0),
phi(getNumberOfArguments(),1.0e-4),
seed(getNumberOfArguments(),time(0))
{
  // Note : parseVector will check that number of arguments is correct
  parseVector("TO",to);
  parseVector("KAPPA",kappa);
  parseVector("PHI",phi);
  parseVector("SEED",seed);
  checkRead();

//  log.printf("  curr");
//  for(unsigned i=0;i<curr.size();i++) log.printf(" %f",curr[i]);
//  log.printf("\n");
//  log.printf("  to");
//  for(unsigned i=0;i<to.size();i++) log.printf(" %f",to[i]);
//  log.printf("\n");
//  log.printf("  with force constant");
//  for(unsigned i=0;i<kappa.size();i++) log.printf(" %f",kappa[i]);
//  log.printf("\n");

  for(unsigned i=0;i<getNumberOfArguments();i++) 
  {
    std::string str_curr=getPntrToArgument(i)->getName()+"_curr";
    addComponent(str_curr); componentIsNotPeriodic(str_curr);
    if(curr[i]!=-1.0) getPntrToComponent(str_curr)->set(curr[i]);
  }
  addComponent("accepted"); componentIsNotPeriodic("accepted");
}

void DIMS::calculate()
{
  bool accepted=true;
  double cv;
  double cv2;
  double k;
  double f;
  float prob;

  for(unsigned i=0;i<getNumberOfArguments();++i)
  {
    cv = difference(i,to[i],getArgument(i));
    cv2=cv*cv;
    k=kappa[i];

    // curr < 0 indicates we have yet to collect a value for it, so we do this
    // now cv2 < curr means that the collective variable is nearer to the
    // target value than it was just previously
    if(curr[i]<0.||cv2<curr[i]) 
    { 
      curr[i] = cv2; 
      accepted = accepted && true
    } 
    // we allow the collective variable to backtrack only with some probability
    else 
    {
      // the probability of backtracking scales inverse exponentially with the
      // amount of backtracking
      prob = std::exp(-1 * std::pow((cv2 - curr[i])/phi[i], 2));
      std::random_device generator(seed[i]);
      std::bernoulli_distribution dist(prob);
      
      // if the backtrack is accepted...
      if dist(generator)
      {
        curr[i] = cv2; 
        accepted = accepted && true;
      }
      // if the backtrack is rejected, apply a harmonic restraining force to
      // the collective variable
      else
      {
        // applying chain rule for the gradient of the potential in CV-space
        // gives us our force on the CV
        f = -2.*k*(cv2-curr[i])*cv;
        setOutputForce(i,f);
        accepted = accepted && false;
      }
    }
    std::string str_curr=getPntrToArgument(i)->getName()+"_curr";
    getPntrToComponent(str_curr)->set(curr[i]);
  }
  getPntrToComponent("accepted")->set(accepted);
}

}
}


