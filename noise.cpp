//by Lin Yangchen

//generates time series of noise of specified length, colour and Gaussian variance
//via spectral synthesis (Cohen et al. 1998 Proc. R. Soc. B 265:11)
//and spectral mimicry (Cohen et al. 1999 Circuits Syst. Signal Process. 18:431).


#include <algorithm>	//std::sort
#include <vector>
#include <fstream>
#include <sys/time.h>

#include <limits>
typedef std::numeric_limits <double> dbl;

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

using namespace std;


/*
==================================================
	SETTINGS
==================================================
*/

//length of time series
const int steps(1000);

//spectral exponent (noise colour)
//white: 0
//pink: 1
//red: 2
//caution: exponents > 2 will produce unrealistic noise
//(increasingly resembling a sine wave)
//using this method of noise generation and are not recommended.
const double fpow(1);

//noise variance
//value <= 0 will generate a vector of zeros
const double env(sqrt(0.05));


/*
==================================================
	PSEUDORANDOM SEED STATE
==================================================
*/

struct timeval uhr;
int systime()
{
	gettimeofday(&uhr, NULL);
	return uhr.tv_usec;
}
const unsigned int seed = systime();
boost::mt19937 rng(seed);	//seed the Mersenne Twister
//boost::mt19937 rng(364527);	//manual seed set


/*
==================================================
	RANDOM NUMBER GENERATORS
==================================================
*/

//uniform random numbers
double runif(double min, double max)
{
	boost::uniform_real <double> u(min, max);
	boost::variate_generator <boost::mt19937&, boost::uniform_real <double> > gen(rng, u);
	return gen();
}

//Gaussian random numbers
double rnorm(double mu, double sd)
{
	boost::variate_generator <boost::mt19937&, boost::normal_distribution <double> >
	gen(rng, boost::normal_distribution <double> (mu, sd));
	return gen();
}


/*
==================================================
	NOISE GENERATOR
==================================================
*/

int main ()
{
	//loop counters
	int timecounter, subtimecounter;


	vector<double> white(steps);
	vector<double> color(steps);
	vector<double> csort(steps);
	vector<double> noise(steps);
	vector<double>::iterator it;


	if(env > 0)
	{
		//define pi
		const double pi(atan(1)*4);

		cout << "noise simulation ..." << endl;

		//phase shifts of sine curves (Cohen et al. 1998)
		double phase[steps/2];
		for(subtimecounter = 0; subtimecounter < steps/2; subtimecounter ++)
		{
			phase[subtimecounter] = runif(0, 2*pi);
		}

		for(timecounter = 0; timecounter < steps; timecounter ++)
		{
			white[timecounter] = rnorm(0, env);

			//spectral synthesis (Cohen et al. 1998)
			color[timecounter] = 0;
			for(subtimecounter = 0; subtimecounter < steps/2; subtimecounter ++)
			{
				color[timecounter] = color[timecounter] +
				sqrt(pow(steps/(2*pi*(subtimecounter + 1)), fpow)) *
				sin(2*pi*(subtimecounter + 1)*(timecounter + 1)/steps +
				    phase[subtimecounter]);
			}
		}

		//spectral mimicry (Cohen et al. 1999):
		//permutate a Gaussian-distributed random sequence of mean 0 and the
		//specified variance such that it matches the spectral density of the above synthesis.
		std::sort(white.begin(), white.end());
		csort = color;
		std::sort(csort.begin(), csort.end());

		for(timecounter = 0; timecounter < steps; timecounter ++)
		{
			it = find(color.begin(), color.end(), csort[timecounter]);
			noise[distance(color.begin(), it)] = white[timecounter];
		}
	} else
	{
		for(timecounter = 0; timecounter < steps; timecounter ++)
		{
			noise[timecounter] = 0;
		}
	}


	//write time series to file
	ofstream fout("timeseries.txt");
	fout.precision(dbl::digits10);	//max precision
	
	for(timecounter = 0; timecounter < steps; timecounter ++)
	{
		fout << scientific << noise[timecounter] << "\n";
	}	
	fout.close();


	//append seed state to file
	ofstream seedout("seed.txt", ios_base::app);
	seedout << seed << endl;
	seedout.close();

	return 0;
}
