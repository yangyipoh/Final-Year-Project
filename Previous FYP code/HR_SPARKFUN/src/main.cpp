//These codes were written to run the whole system

// #include "opencv2/opencv.hpp"
#include "Header.h"

   using namespace std;

   MAX30102 m;

  std::atomic<bool> isQuitRequested;
  std::ofstream fileAcc;
  std::ofstream fileWheel;
  std::ofstream fileBrake;
  std::ofstream fileObject_time;

  std::ofstream rawheartrate;
  std::ofstream rawheartrate_time;
  
  string directory = "/home/pi/Desktop/fyp2022/HR_SPARKFUN/build/Ben/";
  int k = 0;

const string currentTime()
{
	time_t now = time(0);
	struct tm tstruct;
	char buf[80];
	tstruct = *localtime(&now);
	strftime(buf, sizeof(buf), "%Y-%m-%d %X", &tstruct);
	for(char *p = buf; *p != '\0'; p++)
	{
		if(*p == ':')
			*p = ';';
	}
	return buf;
}

 /** @function main */
 int main( int argc, const char** argv )
 {
  //server=init_server();
   cout << "System starting!"  << endl;

	
  wiringPiSetup();

//Wpi12: brake
//Wpi13: acc
//Wpi14: wheel
  
  pinMode(12, OUTPUT);
  pinMode(4, OUTPUT);
  pinMode(14, OUTPUT);

  // Open CV
  pinMode(0, OUTPUT);
  pinMode(1, OUTPUT);

  digitalWrite(1,LOW);
  digitalWrite(0,LOW);
  delay(100);

   // Initialise HR sensor
   gpioTerminate();
   gpioInitialise();
   m.begin();
   m.setup();


   // initialise threads and set them as joinable
   pthread_t hr_coll_thread, mpu;// hr_proc_thread;
   pthread_attr_t attr;
   pthread_attr_init(&attr);
   pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

   // initialise hr data structures
   struct hrthread_data hr_data_coll;
   // arma::vec hr_init_vec = arma::zeros<arma::vec>(30000);

   struct mpu6050 sensor_norm;
   sensor_norm.Norm_brake = Initialize_brake();
   sensor_norm.Norm_acc = Initialize_acc();

// reset isQuitRequested bool (gets set to true when hr data collection is finished)
   isQuitRequested = false;
  

// re-initialise hr data collection struct
   hr_data_coll.hr_reliability = 0;
   hr_data_coll.hr_drowsy = false;
   
   string timeFile = currentTime();
   
      for (int i = 0; i < 5 ; i++)
   {
        digitalWrite(1,HIGH);
  	digitalWrite(0,HIGH);
  	delay(500);  
	digitalWrite(1,LOW);
  	digitalWrite(0,LOW);
  	delay(500);	
   }

   uint32_t start = gpioTick();   
   hr_data_coll.startTime = start;
   sensor_norm.startTime = start;
	
   std::string fileAccName = directory + "MPU6050Acc " + "Name " + "Status" + ".txt";
   std::string fileWheelName = directory + "MPU6050Wheel " + "Name " + "Status" + ".txt";
   std::string fileBrakeName = directory + "MPU6050Brake " + "Name " + "Status" + ".txt";
   std::string fileObject_timeName = directory + "MPU6050output_time "  + "Name " + "Status" + ".txt";

/*
   std::string fileAccName = directory + "MPU6050Acc " + timeFile + ".txt";
   std::string fileWheelName = directory + "MPU6050Wheel " + timeFile + ".txt";
   std::string fileBrakeName = directory + "MPU6050Brake " + timeFile + ".txt";
   std::string fileObject_timeName = directory + "MPU6050output_time "  + timeFile + ".txt";
*/

   fileAcc.open(fileAccName, std::ofstream::out);
   fileWheel.open(fileWheelName, std::ofstream::out);
   fileBrake.open(fileBrakeName , std::ofstream::out);
   fileObject_time.open(fileObject_timeName, std::ofstream::out);
   /*
   rawheartrate.open(directory +  "RawHeartRate " + timeFile + ".txt" , std::ofstream::out);
   rawheartrate_time.open(directory + "RawHeartRate_Time " + timeFile + ".txt" , std::ofstream::out);
   */
   rawheartrate.open(directory +  "RawHeartRate " + "Name " + "Status" + ".txt" , std::ofstream::out);
   rawheartrate_time.open(directory + "RawHeartRate_Time " + "Name " + "Status" + ".txt" , std::ofstream::out);

   pthread_create(&mpu, NULL, MPUthread , &sensor_norm);

// start hr data collection thread
   pthread_create(&hr_coll_thread, NULL, hr_data_collection, (void *)&hr_data_coll);


   void* ret_hr_coll;
   pthread_join(mpu, NULL);
   pthread_join(hr_coll_thread, &ret_hr_coll);

   pthread_attr_destroy(&attr);
   gpioTerminate();
   return 0;
 }


