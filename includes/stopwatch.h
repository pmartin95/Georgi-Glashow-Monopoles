#ifndef _STOPWATCH_H_
#define _STOPWATCH_H_
class mytime{
public:
  void stopwatchStart();
  double stopwatchReadSeconds();
private:
  struct timeval myStartTime;
};
#endif
