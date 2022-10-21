int temp_sensor_1 = A0;
//int temp_sensor_2 = A1;// potentiometer wiper (middle terminal) connected to analog pin 3
                    // outside leads to ground and +5V
int val1 = 0;  // variable to store the value read
//int val2 = 0;

void setup() {
  Serial.begin(9600);           //  setup serial
  pinMode(temp_sensor_1,INPUT);
//  pinMode(temp_sensor_2,INPUT);
}

void loop() {
  val1 = analogRead(temp_sensor_1);  // read the input pin
//  val2 = analogRead(temp_sensor_2);
  Serial.println(val1);          // debug value
  delay(0.1);
}
