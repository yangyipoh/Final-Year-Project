/*********************************************************
Author: Elias Santistevan @ SparkFun Electronics
Date: 5/2021
Library repository can be found here:
https://github.com/sparkfun/SparkFun_KX13X_Arduino_Library

Basic example for reading back accelerometer values using I2C. 

This code is released under the [MIT License](http://opensource.org/licenses/MIT).

Please review the LICENSE.md file included with this example. If you have any questions 
or concerns with licensing, please contact techsupport@sparkfun.com.

Distributed as-is; no warranty is given.
*******************************************************/

#include <Wire.h>
#include "SparkFun_Qwiic_KX13X.h"

//QwiicKX132 kxAccel;
QwiicKX134 kxAccel; // Uncomment this if using the KX134 - check your board
                      //if unsure.  
outputData myData; // This will hold the accelerometer's output. 

float alpha = 0.8;
float gravity[3] = {0, 0, 0}; // gravity board acceleration

// display variables
float acceleration = 0;
float running_avg = 100;

void setup() {
  // Start Serial
  while(!Serial){
    delay(50);
  }
  Serial.begin(9600);
//  Serial.println("Welcome.");
  Wire.begin();
  if( !kxAccel.begin() ){
    Serial.println("Could not communicate with the the KX13X. Freezing.");
    while(1);
  }
//  else
//    Serial.println("Ready.");

  if( !kxAccel.initialize(DEFAULT_SETTINGS)){ // Loading default settings.
    Serial.println("Could not initialize the chip.");
    while(1);
  }
//  else
//    Serial.println("Initialized...");

  // kxAccel.setRange(KX132_RANGE16G);
  // kxAccel.setRange(KX134_RANGE32G); // For a larger range uncomment

}

void loop() {
  acceleration = 0;
  
  for (int i=0; i<running_avg; ++i) {
    myData = kxAccel.getAccelData(); 

    float x = myData.xData;
    float y = myData.yData;
    float z = myData.zData;
    
    gravity[0] = alpha*gravity[0] + (1-alpha)*x;
    gravity[1] = alpha*gravity[1] + (1-alpha)*y;
    gravity[2] = alpha*gravity[2] + (1-alpha)*z;
  
    x = x - gravity[0];
    y = y - gravity[1];
    z = z - gravity[2];
  
    acceleration += sqrt(x*x + y*y + z*z);
  }

  acceleration /= running_avg;
//  Serial.print("X: ");
//  Serial.print(x, 4);
//  Serial.print("g ");
//  Serial.print(" Y: ");
//  Serial.print(y, 4);
//  Serial.print("g ");
//  Serial.print(" Z: ");
//  Serial.print(z, 4);
//  Serial.println("g ");
  Serial.print("Acceleration: ");
  Serial.print(acceleration*100, 4);
  Serial.print('\n');

  delay(20); 

}
