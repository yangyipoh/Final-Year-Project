//These codes were written for the operation of MPU6050 sensors
//to collect angle data

#include "Header.h"

//
//MPU6050_startpart1
//

#include <wiringPiI2C.h>
#include <math.h>

#define PWR_MGMT_1   0x6B
#define SMPLRT_DIV   0x19
#define CONFIG       0x1A
#define GYRO_CONFIG  0x1B
#define INT_ENABLE   0x38

using namespace std;
using namespace std::chrono;

milliseconds elapsed; //time elapsed

int first_sr_cap, last_sr_cap;
int fd;
double value;
int acclX, acclY, acclZ;
double acclX_scaled, acclY_scaled, acclZ_scaled;
int sensor,i;

extern int k;
extern std::atomic<bool> isQuitRequested;

  //Wpi12: brake
  //Wpi13: acceleration
  //Wpi14: wheel

void MPU6050_Init(){
	
wiringPiI2CWriteReg8 (fd, SMPLRT_DIV, 0x07);	/* Write to sample rate register */
wiringPiI2CWriteReg8 (fd, PWR_MGMT_1, 0x01);	/* Write to power management register */
wiringPiI2CWriteReg8 (fd, CONFIG, 0);		/* Write to Configuration register */
wiringPiI2CWriteReg8 (fd, GYRO_CONFIG, 24);	/* Write to Gyro Configuration register */
wiringPiI2CWriteReg8 (fd, INT_ENABLE, 0x01);	/*Write to interrupt enable register */


} 

  //
  //Calculation of x, y, z axis rotation angle of sensor  
  //

int read_word_2c(int addr, int fd)
{
  int val;
  val = wiringPiI2CReadReg8(fd, addr);
  val = val << 8;
  val += wiringPiI2CReadReg8(fd, addr+1);
  if (val >= 0x8000)
    val = -(65536 - val);

  return val;
}

inline double dist(double a, double b)
{
  return sqrt((a*a) + (b*b));
}

inline double get_y_rotation(double x, double y, double z)
{
  double radians;
  radians =1.13* atan2(y, dist(x, z));
  return -(radians * (180.0 / M_PI));
}

inline double get_x_rotation(double x, double y, double z)
{
  double radians;
  radians = 1.13* atan2(x, dist(y, z));
  return (radians * (180.0 / M_PI));
}

inline double get_z_rotation(double x, double y, double z)
{
  double radians;
  radians = 1.13* atan2(dist(x,y), z);
  return (radians * (180.0 / M_PI));
}

  //Function to calculate initial brake pedal angle
double Initialize_brake(void)
{
  double Sum = 0;
  digitalWrite(12,HIGH);
  delay(10);
  fd = wiringPiI2CSetup (0x69);
  MPU6050_Init();
  for(i=0;i<50;i++)
  {
    acclX = read_word_2c(0x3B,fd);
    acclY = read_word_2c(0x3D,fd);
    acclZ = read_word_2c(0x3F,fd);

    acclX_scaled = acclX / 16384.0;
    acclY_scaled = acclY / 16384.0;
    acclZ_scaled = acclZ / 16384.0;
    Sum += get_y_rotation(acclX_scaled, acclY_scaled, acclZ_scaled);
    delay(20);
  }
  digitalWrite(12,LOW);
  printf("\nInitial brake value is: %.2f\n",Sum/50);
  return Sum/50;
}

  //Function to calculate initial acceleration pedal angle
double Initialize_acc(void)
{
  double Sum = 0;
  digitalWrite(4,HIGH);
  delay(10);
  fd = wiringPiI2CSetup (0x69);
  MPU6050_Init();
  for(i=0;i<50;i++)
  {
    acclX = read_word_2c(0x3B,fd);
    acclY = read_word_2c(0x3D,fd);
    acclZ = read_word_2c(0x3F,fd);

    acclX_scaled = acclX / 16384.0;
    acclY_scaled = acclY / 16384.0;
    acclZ_scaled = acclZ / 16384.0;
    Sum += get_y_rotation(acclX_scaled, acclY_scaled, acclZ_scaled);
  }
  digitalWrite(4,LOW);
  printf("Initial accelerator value is: %.2f\n\n",Sum/50);
  return Sum/50;
}

  //Create text files to store angle data along with time of 
  //accelearation pedal, brake pedal and steering wheel
  extern ofstream fileAcc;
  extern ofstream fileWheel;
  extern ofstream fileBrake;
  extern ofstream fileObject_time;
  double brakePedal, accPedal, wheel;

//
//MPU6050_endpart1
//

void *MPUthread(void* Input)
{
	struct mpu6050 *sensorptr = (struct mpu6050 *) Input;
	
	//
	//MPU6050_startpart2
	//	
	
	fileObject_time << gpioTick() - sensorptr->startTime << endl; //Record starting time

	while ( isQuitRequested == false)   // time limit
	{
		for (sensor = 1; sensor < 4 ; sensor++)  //Read values from each sensor consecutively
  		{					 //by setting its AD0 pin high and low
        switch (sensor)
   			{
          case 1:{
            digitalWrite(12,HIGH);
            break;
          }
          case 2:{
            digitalWrite(4,HIGH);
            break;
          }  
          case 3:{
            digitalWrite(14,HIGH);
            break;
          }        
	   		}

        fd = wiringPiI2CSetup (0x69);
        MPU6050_Init();
  

        acclX = read_word_2c(0x3B,fd);
        acclY = read_word_2c(0x3D,fd);
        acclZ = read_word_2c(0x3F,fd);

        acclX_scaled = acclX / 16384.0;
        acclY_scaled = acclY / 16384.0;
        acclZ_scaled = acclZ / 16384.0;
    
         
        if(sensor == 1)	//Record values from the sensor attached to brake pedal
        {
          brakePedal = round(100.0 * (get_y_rotation(acclX_scaled, acclY_scaled, acclZ_scaled) - sensorptr->Norm_brake)) / 100.0;
          fileBrake << brakePedal << std::endl;
          //std::cout << "BRAKE: " << brakePedal << std::endl;
          delay(10);
        }

        if(sensor == 2)  //Record values from the sensor attached to acceleration pedal
        {
	
          accPedal = round(100.0 * (get_y_rotation(acclX_scaled, acclY_scaled, acclZ_scaled) - sensorptr->Norm_acc)) / 100.0;
          fileAcc << accPedal << std::endl;
          // std::cout << "ACC: " << accPedal <<std::endl;
          delay(10);
        }

        if(sensor == 3)  //Record values from the sensor attached to steerign wheel
        {
          wheel = round(100.0 * get_y_rotation(acclX_scaled, acclY_scaled, acclZ_scaled)) / 100.0;
          fileWheel << wheel << std::endl;
          std::cout << "Wheel: " << wheel <<std::endl;
          fileObject_time << gpioTick() - sensorptr->startTime << std::endl; //Record the time of this one cycle
          delay(10);
        }

   
        switch (sensor)
   			{
          case 1:{
            digitalWrite(12,LOW);
            break;
          }
          case 2:{
            digitalWrite(4,LOW);
            break;
          }  
          case 3:{
            digitalWrite(14,LOW);
            break;
          }        
   			}//endswitch
  		}//endforloop
	}//endwhileloop
	fileObject_time  << gpioTick() - sensorptr->startTime << endl; //Record the ending time
	fileAcc.close();
  fileWheel.close();
  fileBrake.close();
  fileObject_time.close();
  pthread_exit(NULL);
}//end thread function

//
//MPU6050_endpart2
//
