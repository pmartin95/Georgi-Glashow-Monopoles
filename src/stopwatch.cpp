#include <cstdlib>
#include <sys/time.h>
#include "stopwatch.h"



void mytime::stopwatchStart()
{
								gettimeofday(&myStartTime, NULL);
}

double mytime::stopwatchReadSeconds()
{
								struct timeval endTime;
								gettimeofday(&endTime, 0);

								long ds = endTime.tv_sec - myStartTime.tv_sec;
								long dus = endTime.tv_usec - myStartTime.tv_usec;
								return ds + 0.000001*dus;
}
