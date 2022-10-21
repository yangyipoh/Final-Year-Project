// libraries
#include <Wire.h>
#include "SparkFun_Qwiic_KX13X.h"

// variables
int acc_main = 13;
int acc_ref = 12;
QwiicKX134 kxAccel_main;
QwiicKX134 kxAccel_ref;  
outputData myData_main;
outputData myData_ref;

// ------------------------ ISR setup -----------------------------------
void isr_setup() {
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

// ------------------------- Read Accelerometer Data --------------------
outputData read_acc_data(QwiicKX134* kxAccel) {
  return kxAccel->getAccelData();
}

void print_data(outputData* data_main, outputData* myData_ref) {
  Serial.print(data_main->xData, 4);
  Serial.print(',');
  Serial.print(data_main->yData, 4);
  Serial.print(',');
  Serial.print(data_main->zData, 4);
  Serial.print(',');

  Serial.print(myData_ref->xData, 4);
  Serial.print(',');
  Serial.print(myData_ref->yData, 4);
  Serial.print(',');
  Serial.println(myData_ref->zData, 4);
}

// ------------------------ Accelerometer setup -------------------------

void acc_setup(QwiicKX134* kxAccel, int idx) {
  // Main accelerometer --> default address
  if (idx == 1){
    if( !kxAccel->begin() ){
      Serial.println("Could not communicate with the the KX134_Main. Freezing.");
      while(1);
    }
  }
  // Ref accelerometer --> Alternative address
  else {
    if( !kxAccel->begin(KX13X_ALT_ADDRESS) ){
      Serial.println("Could not communicate with the the KX134_Ref. Freezing.");
      while(1);
    }
  }
  
  if( !kxAccel->initialize(DEFAULT_SETTINGS)){ // Loading default settings.
    Serial.println("Could not initialize the chip.");
    while(1);
  }
  // kxAccel.setRange(KX132_RANGE16G);
  // kxAccel.setRange(KX134_RANGE32G); // For a larger range uncomment
}

void setup(){
  // start serial
  while(!Serial){
    delay(50);
  }
  Serial.begin(115200);

  Wire.begin();
  
  // pin setup
  pinMode(acc_main, OUTPUT);
  pinMode(acc_ref, OUTPUT);
  
  // setup accelerometer
  acc_setup(&kxAccel_main, 1);
  acc_setup(&kxAccel_ref, 2);

  // setup timer interrupt @ 50Hz 
  isr_setup();
}

// --------------------- ISR -------------------------------
ISR(TIMER4_COMPA_vect){
  sei();
  // main acceleration data
  myData_main = read_acc_data(&kxAccel_main);
  myData_ref = read_acc_data(&kxAccel_ref);
  
  print_data(&myData_main, &myData_ref);
}
 
void loop(){
  //do other things here
  delay(5000);
}
