// MC simulation for Pixi detector
// Batch version
//==============================================================


#include "MC.h"
#include "TExperiment.h"

TRandom3 *rnd;
TExperiment *Experiment;

int main(int argn, char *argv[])
{
  // open option file and get all parameters
  /* this should go in a dedicated class */
  if (argn != 2){
    std::cout<< "Error: a config file expected! Exiting..."<< std::endl;
    return 1;
  }
  std::ifstream config_file(argv[1], ios::in);
  if( !config_file )
    {
      std::cout<< "Error: config file cannot be opened! Exiting..."<< std::endl;
      return 1;
    }
  std::cout<< "Parsing config file "<< argv[1] << std::endl;
  // parse the file
  std::string sLine, sKeyWord, sValue;
  const string  sDelim( ":" );
  std::size_t delim_pos = 0 ;
  // Map to associate the strings with the enum values
  std::map<std::string, int> s_mapStringValues;
  s_mapStringValues["gas_id"]   = 1;
  s_mapStringValues["pressure"] = 2;
  s_mapStringValues["thickness"] = 3;
  s_mapStringValues["num_evts"] = 4;
  s_mapStringValues["energy_min"] = 5;
  s_mapStringValues["energy_max"] = 6;
  s_mapStringValues["pol_angle"] = 7;
  s_mapStringValues["pol_degree"] = 8;
  int gas_id = -1;
  double pressure = -1.;
  double thickness = -1.;
  int num_evts = -1;
  double energy_min = -1.;
  double energy_max = -1.;
  double pol_angle = -1.;
  double pol_degree = -1.;
  while(std::getline(config_file,sLine))
    {
      if(!sLine.empty() && sLine.find("#")!=0) // empy line or start with #
	{
	  delim_pos = sLine.find_first_of( sDelim );// find delim
	  sKeyWord  = sLine.substr( 0, delim_pos ); // extract key word
	  sValue    = sLine.substr(delim_pos+1);    // extract value
	  
	  //std::cout << "Line(" << sLine << ") " << sKeyWord
	  //<< "-" << sValue<< std::endl;
	  switch (s_mapStringValues[sKeyWord])
	    {
	    case 1:
	      gas_id = atoi( sValue.c_str());
	      std::cout <<"Gas Id set to "<< gas_id <<std::endl;
	      break;
	    case 2:
	      //pressure = std::stod(sValue); // requires c++11
	      pressure = atof( sValue.c_str());
	      std::cout <<"Pressure set to "<< pressure <<std::endl;
	      break;
	    case 3:
	      thickness = atof( sValue.c_str());
	      std::cout <<"thickness set to "<< thickness <<std::endl;
	      break;
	    case 4:
	      num_evts = atoi( sValue.c_str());
	      std::cout <<"num_evts set to "<< num_evts <<std::endl;
	      break;
	    case 5:
	      energy_min = atof( sValue.c_str());
	      std::cout <<"energy_min set to "<< energy_min <<std::endl;
	      break;
	    case 6:
	      energy_max = atof( sValue.c_str());
	      std::cout <<"energy_max set to "<< energy_max <<std::endl;
	      break;
	    case 7:
	      pol_angle = atof( sValue.c_str());
	      std::cout <<"pol_angle set to "<< pol_angle  <<std::endl;
	      break;
	    case 8:
	      pol_degree = atof( sValue.c_str());
	      std::cout <<"pol_degree set to "<< pol_degree <<std::endl;
	      break;
	    default:
	      std::cout << "value of "<< sKeyWord <<" unknown"<<std::endl;
	    }
	}
    }
 
  //
  // running the simulation...
  //
  std::cout<< "Starting MC in batch mode" << std::endl;

  
  rnd = new TRandom3();
  Experiment = new TExperiment(rnd);
  Experiment->SetMixID(gas_id);
  Experiment->SetPressure(pressure);
  Experiment->SetThickness(thickness);
  Experiment->Generalsetup();
  // source TBD, now is 'monochromatic' by default
  // Using only Energy Min as energy. scan tbd!!!
  Experiment->SetSource(pol_angle,pol_degree, energy_min);
  Experiment->EventsTree(num_evts);
  
  std::cout<< "End of run!" << std::endl;
  return 0;
}

