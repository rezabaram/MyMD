// This file is a part of Molecular Dynamics code for
// simulating ellipsoidal packing. The author cannot
// guarantee the correctness nor the intended functionality.
//
// March 2012, Reza Baram

#include"include/main.h"
#include"include/version.h"
#include<fstream>
#include<ctime>
#include<cstring>

using namespace std;


long RNGSeed;
extern MTRand rgen;

string config_file="config";
CConfig config;

// The seed is scrambled before use.  This is the original formula, kept
// exactly: changing it would change every result in the regression references
// without changing any physics.
static long scramble_seed(long seed){ return 313*seed + 1; }

static void print_usage(ostream &out, const char *prog)
{
	out<<
"Usage: "<<prog<<" [options] [seed] [config-file]\n"
"\n"
"Runs a simulation.  The two positional arguments are kept for compatibility\n"
"with the original '"<<prog<<" <seed> <config-file>' form, which bin/run.sh\n"
"and the Makefile use.\n"
"\n"
"Options:\n"
"  -c, --config FILE     configuration file (default: config)\n"
"  -s, --seed N          random seed (default: 0)\n"
"  -D, --set KEY=VALUE   override a parameter; may be repeated\n"
"  -o, --output PREFIX   override the 'output' parameter\n"
"      --print-config    print the effective parameters and exit\n"
"      --save-config F   write the effective parameters to F\n"
"                        (default: config.used in the working directory)\n"
"      --no-save-config  do not write config.used\n"
"  -h, --help            this text\n"
"  -V, --version         version\n"
"\n"
"Examples:\n"
"  "<<prog<<" 3 config_quick\n"
"  "<<prog<<" -c config_stillinger -s 7 -D nParticle=500 -D maxTime=2\n"
"  "<<prog<<" --print-config -c config_quick\n"
	;
}

// Writes the parameters actually in force, with a header recording where they
// came from, so a result can be traced back to the inputs that produced it.
static void save_effective_config(const string &path)
{
	ofstream out(path.c_str());
	if(!out.good()){
		WARNING("could not write "<<path<<"; continuing without it");
		return;
		}
	time_t now=time(NULL);
	char stamp[64];
	struct tm *lt=localtime(&now);
	if(!lt || strftime(stamp, sizeof(stamp), "%Y-%m-%dT%H:%M:%S", lt)==0)
		stamp[0]='\0';
	out<<"# ellipmd "<<MYMD_VERSION<<" effective configuration\n"
	   <<"# seed "<<RNGSeed<<"\n"
	   <<"# config file "<<config_file<<"\n"
	   <<"# written "<<stamp<<"\n";
	config.print(out);
	}


void Initialize(){
	rgen.seed(RNGSeed);
	config.parse(config_file);
	}

void Run(){

	CSys sys(config.get_param<unsigned int>("nParticle"));
	sys.initialize(config);
	sys.solve();
	}

int main(int argc, char **params){
	string save_config="config.used";
	bool write_config=true;
	bool print_only=false;
	string output_override;
	vector<string> overrides;

	// ---- parse arguments -------------------------------------------------
	// Positional arguments are still accepted, so bin/run.sh and anything else
	// written against the original interface keeps working.
	vector<string> positional;
	for(int i=1; i<argc; ++i){
		string a=params[i];
		string value;
		bool has_value=false;
		if(a.size()>2 && a[0]=='-' && a[1]=='-'){
			string::size_type eq=a.find('=');
			if(eq!=string::npos){
				value=a.substr(eq+1);
				a=a.substr(0,eq);
				has_value=true;
				}
			}
		#define NEED_VALUE(opt) \
			if(!has_value){ \
				if(i+1>=argc){ \
					cerr<<(opt)<<" needs an argument"<<endl; \
					return 2; \
					} \
				value=params[++i]; \
				}

		if(a=="-h" || a=="--help"){ print_usage(cout, params[0]); return 0; }
		else if(a=="-V" || a=="--version"){
			cout<<"ellipmd "<<MYMD_VERSION<<endl;
			return 0;
			}
		else if(a=="-c" || a=="--config"){ NEED_VALUE(a); config_file=value; }
		else if(a=="-s" || a=="--seed"){
			NEED_VALUE(a);
			RNGSeed=scramble_seed(atol(value.c_str()));
			}
		else if(a=="-o" || a=="--output"){ NEED_VALUE(a); output_override=value; }
		else if(a=="-D" || a=="--set"){ NEED_VALUE(a); overrides.push_back(value); }
		else if(a=="--save-config"){ NEED_VALUE(a); save_config=value; }
		else if(a=="--no-save-config"){ write_config=false; }
		else if(a=="--print-config"){ print_only=true; }
		else if(a.size()>1 && a[0]=='-'){
			cerr<<"unknown option: "<<a<<"\n\n";
			print_usage(cerr, params[0]);
			return 2;
			}
		else positional.push_back(a);
		#undef NEED_VALUE
		}

	// positional seed and config file, in the original order
	if(positional.size()>2){
		cerr<<"too many arguments; expected at most a seed and a config file"<<endl;
		return 2;
		}
	if(positional.size()>=1) RNGSeed=scramble_seed(atol(positional[0].c_str()));
	if(positional.size()>=2) config_file=positional[1];

	try {
		rgen.seed(RNGSeed);
		config.parse(config_file);

		for(size_t i=0; i<overrides.size(); ++i){
			string::size_type eq=overrides[i].find('=');
			if(eq==string::npos){
				cerr<<"--set expects KEY=VALUE, got '"<<overrides[i]<<"'"<<endl;
				return 2;
				}
			config.set(overrides[i].substr(0,eq), overrides[i].substr(eq+1));
			}
		if(!output_override.empty())
			config.set("output", output_override);

		if(print_only){
			cout<<"# ellipmd "<<MYMD_VERSION
			    <<"   seed "<<RNGSeed
			    <<"   config "<<config_file<<endl;
			config.print(cout);
			return 0;
			}

		if(write_config) save_effective_config(save_config);

		cerr<<"ellipmd "<<MYMD_VERSION<<", RNG seed "<<RNGSeed
		    <<", config "<<config_file<<endl;
		Run();
		return 0;
	} catch(CException &e)
	{
	e.Report();
	return 1;
	}
}
