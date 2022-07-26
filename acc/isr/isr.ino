// definition
int temp_sensor_1 = A0;
int val1 = 0;

// libraries
#include <Wire.h>
#include "SparkFun_Qwiic_KX13X.h"

// variables
QwiicKX134 kxAccel;  
outputData myData;

//float alpha = 0.8;
//float gravity[3] = {0, 0, 0}; // gravity board acceleration

// display variables
//float acceleration = 0;

void setup(){
  // ----------------- setting up accelerometer --------------------------
  while(!Serial){
    delay(50);
  }
  Serial.begin(115200);
  Wire.begin();
  if( !kxAccel.begin() ){
    Serial.println("Could not communicate with the the KX13X. Freezing.");
    while(1);
  }

  if( !kxAccel.initialize(DEFAULT_SETTINGS)){ // Loading default settings.
    Serial.println("Could not initialize the chip.");
    while(1);
  }
  pinMode(temp_sensor_1,INPUT);

  // kxAccel.setRange(KX132_RANGE16G);
  // kxAccel.setRange(KX134_RANGE32G); // For a larger range uncomment

  // ------------------ setting up timer interrupt ------------------------- 
  cli();//stop interrupts
  //set timer4 interrupt at 1Hz
  TCCR4A = 0;// set entire TCCR1A register to 0
  TCCR4B = 0;// same for TCCR1B
  TCNT4  = 0;//initialize counter value to 0
  // set compare match register for 50hz increments
  OCR4A = 1249;// = (16*10^6) / (256*50) - 1 (must be <65536)
  // turn on CTC mode
  TCCR4B |= (1 << WGM12);
  // Set CS12 bits for 256 prescaler
  TCCR4B |= (1 << CS12);
  // enable timer compare interrupt
  TIMSK4 |= (1 << OCIE4A);
  sei();//allow interrupts
}

// ------------------- ISR -------------------------------
ISR(TIMER4_COMPA_vect){
  sei();
  // acceleration data
  myData = kxAccel.getAccelData(); 

  float x = myData.xData;
  float y = myData.yData;
  float z = myData.zData;
  
  Serial.print(x, 4);
  Serial.print(',');
  Serial.print(y, 4);
  Serial.print(',');
  Serial.print(z, 4);
  Serial.print(',');

  // temperature data
  val1 = analogRead(temp_sensor_1);
  Serial.println(val1);
}
 
void loop(){
  //do other things here
  delay(1000);
}
