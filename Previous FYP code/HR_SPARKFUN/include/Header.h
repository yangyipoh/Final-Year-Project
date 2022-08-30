#ifndef HEADER
#define HEADER


#include <iostream>
#include <stdio.h>
#include <chrono>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstring>

#include <string>
#include <wiringPi.h>

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/socket.h>


#include <pigpio.h>
#include "MAX30105.h"

#include <pthread.h>
#include <atomic>

void head_stats(int angle);
void *hr_data_collection(void* threadarg);
void *MPUthread(void* Input);
void *threadVideoRecord(void *f);
double Initialize_acc(void);
double Initialize_brake(void);
void MPU6050_Init();
const std::string currentTime();

/** Struct defines */
struct hrthread_data {
	//arma::vec hr;
	//arma::vec signal_quality20secs;
	//arma::vec minuteIndexes;
	float hr_reliability;
	bool hr_drowsy;
	float min_count;
	uint32_t startTime;
};

struct mpu6050 {
	double Norm_brake;
	double Norm_acc;
	uint32_t startTime;

};

 extern MAX30102 m;
 extern std::ofstream myfile;
 extern std::ifstream datafile;
 extern std::ofstream hrtime;
 extern std::ofstream framerate;

#endif



