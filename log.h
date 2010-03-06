//http://www.nongnu.org/svas/log_8h-source.html
#ifndef LOG_H
#define LOG_H 
#include <unistd.h>
#include <time.h>
#include <stdio.h>
#include <iostream>
#include <string>

#define log std::cout<<TS<<"  "<<__FILE__<<":"<<__LINE__<<"  "

#define TS timeStamp()

std::string timeStamp(void);

unsigned long int sTime(void);


void andTimePasses(void);

unsigned long int localTime=0;

unsigned long int sTime(void){ 
return localTime; 
};

std::string timeStamp(){
char str[256];
sprintf(str,"%d::%ld",getpid(),sTime());

return std::string(str);
}

void andTimePasses(){
localTime++;
}

#endif /* LOG_H */
