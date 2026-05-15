/*
  ESP32 EKF SOC Estimator with 1-RC Battery Model
  ------------------------------------------------
  Translated from Python EKF simulation.

  State vector:
    x[0] = SOC, range 0.0 to 1.0
    x[1] = VRC, transient RC polarization voltage

  Measurement model:
    V_terminal = OCV(SOC) - I*R0 - VRC

  Notes:
    - Positive current is assumed to mean DISCHARGING.
    - If SOC goes up while discharging, flip CURRENT_SIGN to -1.0.
    - OCV lookup is linear interpolation for Arduino simplicity.
    - You can copy CSV output from Serial Monitor into Excel,
      or use the Serial Plotter mode for live plotting.
*/

#include <Wire.h>
#include <INA226_WE.h>
#include <math.h>

// -------------------- Hardware setup --------------------
#define INA226_ADDRESS 0x40

#define SDA_PIN 21
#define SCL_PIN 22

#define MOSFET_CRITICAL   25
#define MOSFET_MEDIUM     26
#define MOSFET_DEFERRABLE 27

INA226_WE ina226 = INA226_WE(INA226_ADDRESS);

// -------------------- Output mode --------------------
// true  = Arduino Serial Plotter format
// false = CSV format for Excel logging
const bool SERIAL_PLOTTER_MODE = false;

// -------------------- Timing --------------------
const unsigned long EKF_INTERVAL_MS = 500;    // EKF update every 0.5s
const unsigned long LOG_INTERVAL_MS = 5000;   // Excel/CSV log every 5s

unsigned long lastEKFTime = 0;
unsigned long lastLogTime = 0;
unsigned long startTime = 0;

// -------------------- Battery / model parameters --------------------
// Replace these with your measured values later.
const float BATTERY_CAPACITY_AH = 0.820;                 // Ah
const float CAPACITY_AS = BATTERY_CAPACITY_AH * 3600.0;  // Amp-seconds

const float R0 = 0.020;    // Ohms, instantaneous internal resistance
const float R1 = 0.030;    // Ohms, RC branch resistance
const float C1 = 2000.0;   // Farads, RC branch capacitance

// If INA226 current is negative during discharge, use -1.0 here.
const float CURRENT_SIGN = 1.0;

// -------------------- EKF tuning values --------------------
// Initial covariance P
float P00 = 0.01;
float P01 = 0.0;
float P10 = 0.0;
float P11 = 0.001;

// Process noise Q
const float Q00 = 1e-8;
const float Q01 = 0.0;
const float Q10 = 0.0;
const float Q11 = 1e-5;

// Measurement noise R
// Python used R = 0.0064, which means voltage std dev around 0.08 V.
// Tune this based on how noisy your voltage measurement is.
const float R_MEAS = 0.0064;

// -------------------- EKF state --------------------
float soc = 0.90;   // Initial SOC estimate, 0.0 to 1.0
float vrc = 0.0;    // Initial RC transient voltage estimate

// -------------------- OCV table --------------------
// Same values from Python notebook.
const int OCV_N = 21;

const float socPoints[OCV_N] = {
  0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50,
  55, 60, 65, 70, 75, 80, 85, 90, 95, 100
};

const float ocvPoints[OCV_N] = {
  3.00, 3.13, 3.40, 3.55, 3.62, 3.66, 3.68, 3.70,
  3.71, 3.72, 3.73, 3.74, 3.75, 3.77, 3.80, 3.84,
  3.90, 3.98, 4.06, 4.14, 4.20
};

// -------------------- Helper functions --------------------
float clampFloat(float value, float minVal, float maxVal) {
  if (value < minVal) return minVal;
  if (value > maxVal) return maxVal;
  return value;
}


float ocvTangents[OCV_N];
// Computes monotonic cubic Hermite tangents for OCV curve.
// Run this once in setup().
void computeOCVTangents() {
  float slopes[OCV_N - 1];

  // Secant slopes between each pair of points
  for (int i = 0; i < OCV_N - 1; i++) {
    slopes[i] = (ocvPoints[i + 1] - ocvPoints[i]) / (socPoints[i + 1] - socPoints[i]);
  }

  // Endpoint tangents
  ocvTangents[0] = slopes[0];
  ocvTangents[OCV_N - 1] = slopes[OCV_N - 2];

  // Interior tangents
  for (int i = 1; i < OCV_N - 1; i++) {
    if (slopes[i - 1] * slopes[i] <= 0.0) {
      ocvTangents[i] = 0.0;
    } else {
      ocvTangents[i] = 0.5 * (slopes[i - 1] + slopes[i]);
    }
  }

  // Monotonic limiting step to reduce overshoot
  for (int i = 0; i < OCV_N - 1; i++) {
    if (slopes[i] == 0.0) {
      ocvTangents[i] = 0.0;
      ocvTangents[i + 1] = 0.0;
    } else {
      float a = ocvTangents[i] / slopes[i];
      float b = ocvTangents[i + 1] / slopes[i];

      float magnitude = sqrt(a * a + b * b);

      if (magnitude > 3.0) {
        float scale = 3.0 / magnitude;
        ocvTangents[i] = scale * a * slopes[i];
        ocvTangents[i + 1] = scale * b * slopes[i];
      }
    }
  }
}

// Monotonic cubic Hermite spline OCV(SOC).
// Input SOC is 0.0 to 1.0.
// Output voltage in volts.
float getOCV(float soc01) {
  float socPercent = clampFloat(soc01 * 100.0, 0.0, 100.0);

  if (socPercent <= socPoints[0]) return ocvPoints[0];
  if (socPercent >= socPoints[OCV_N - 1]) return ocvPoints[OCV_N - 1];

  for (int i = 0; i < OCV_N - 1; i++) {
    if (socPercent >= socPoints[i] && socPercent <= socPoints[i + 1]) {
      float x0 = socPoints[i];
      float x1 = socPoints[i + 1];
      float y0 = ocvPoints[i];
      float y1 = ocvPoints[i + 1];

      float m0 = ocvTangents[i];
      float m1 = ocvTangents[i + 1];

      float h = x1 - x0;
      float t = (socPercent - x0) / h;

      float t2 = t * t;
      float t3 = t2 * t;

      float h00 = 2.0 * t3 - 3.0 * t2 + 1.0;
      float h10 = t3 - 2.0 * t2 + t;
      float h01 = -2.0 * t3 + 3.0 * t2;
      float h11 = t3 - t2;

      return h00 * y0 + h10 * h * m0 + h01 * y1 + h11 * h * m1;
    }
  }

  return ocvPoints[OCV_N - 1];
}


// Derivative dOCV/dSOC, where SOC is 0.0 to 1.0.
// Since the spline x-axis is SOC percent, multiply by 100 at the end.
float getdOCVdSOC(float soc01) {
  float socPercent = clampFloat(soc01 * 100.0, 0.0, 100.0);

  if (socPercent <= socPoints[0]) {
    return ocvTangents[0] * 100.0;
  }

  if (socPercent >= socPoints[OCV_N - 1]) {
    return ocvTangents[OCV_N - 1] * 100.0;
  }

  for (int i = 0; i < OCV_N - 1; i++) {
    if (socPercent >= socPoints[i] && socPercent <= socPoints[i + 1]) {
      float x0 = socPoints[i];
      float x1 = socPoints[i + 1];
      float y0 = ocvPoints[i];
      float y1 = ocvPoints[i + 1];

      float m0 = ocvTangents[i];
      float m1 = ocvTangents[i + 1];

      float h = x1 - x0;
      float t = (socPercent - x0) / h;

      float t2 = t * t;

      float dh00 = 6.0 * t2 - 6.0 * t;
      float dh10 = 3.0 * t2 - 4.0 * t + 1.0;
      float dh01 = -6.0 * t2 + 6.0 * t;
      float dh11 = 3.0 * t2 - 2.0 * t;

      // This derivative is currently V per SOC percent.
      float dOCV_dPercent =
        (dh00 * y0 + dh10 * h * m0 + dh01 * y1 + dh11 * h * m1) / h;

      // Convert from V per percent to V per SOC unit.
      return dOCV_dPercent * 100.0;
    }
  }

  return 0.0;
}

// Battery terminal voltage prediction.
float predictVoltage(float soc01, float currentA, float vrcEstimate) {
  return getOCV(soc01) - (currentA * R0) - vrcEstimate;
}

// Load shedding logic based on estimated SOC.
// Adjust thresholds however you want.
void updateLoads(float socEstimate) {
  digitalWrite(MOSFET_CRITICAL, HIGH); // critical load always on

  if (socEstimate > 0.70) {
    digitalWrite(MOSFET_MEDIUM, HIGH);
    digitalWrite(MOSFET_DEFERRABLE, HIGH);
  }
  else if (socEstimate > 0.40) {
    digitalWrite(MOSFET_MEDIUM, HIGH);
    digitalWrite(MOSFET_DEFERRABLE, LOW);
  }
  else {
    digitalWrite(MOSFET_MEDIUM, LOW);
    digitalWrite(MOSFET_DEFERRABLE, LOW);
  }
}

void runEKF(float measuredVoltageV, float measuredCurrentA, float dt) {
  // -------------------- Predict step --------------------
  float a22 = exp(-dt / (R1 * C1));

  float socPred = soc - ((measuredCurrentA * dt) / CAPACITY_AS);
  socPred = clampFloat(socPred, 0.0, 1.0);

  float vrcPred = (a22 * vrc) + (R1 * (1.0 - a22) * measuredCurrentA);

  // A = [1, 0; 0, a22]
  // pred_error = A @ P @ A.T + Q
  float Ppred00 = P00 + Q00;
  float Ppred01 = P01 * a22 + Q01;
  float Ppred10 = P10 * a22 + Q10;
  float Ppred11 = (a22 * a22 * P11) + Q11;
  

  // -------------------- Measurement update --------------------
  // H = [dOCV/dSOC, -1]
  float H0 = getdOCVdSOC(socPred);
  float H1 = -1.0;

  float predVoltage = predictVoltage(socPred, measuredCurrentA, vrcPred);
  float innovation = measuredVoltageV - predVoltage;

  // S = H * Ppred * H^T + R
  float S =
    H0 * (Ppred00 * H0 + Ppred01 * H1) +
    H1 * (Ppred10 * H0 + Ppred11 * H1) +
    R_MEAS;

  // Avoid division by zero or unstable tiny S.
  if (fabs(S) < 1e-9) {
    return;
  }

  // K = Ppred * H^T / S
  float K0 = (Ppred00 * H0 + Ppred01 * H1) / S;
  float K1 = (Ppred10 * H0 + Ppred11 * H1) / S;

  // x = xPred + K * innovation
  soc = socPred + K0 * innovation;
  vrc = vrcPred + K1 * innovation;

  soc = clampFloat(soc, 0.0, 1.0);

  // P = (I - K*H) * Ppred
  float IKH00 = 1.0 - K0 * H0;
  float IKH01 = 0.0 - K0 * H1;
  float IKH10 = 0.0 - K1 * H0;
  float IKH11 = 1.0 - K1 * H1;

  float newP00 = IKH00 * Ppred00 + IKH01 * Ppred10;
  float newP01 = IKH00 * Ppred01 + IKH01 * Ppred11;
  float newP10 = IKH10 * Ppred00 + IKH11 * Ppred10;
  float newP11 = IKH10 * Ppred01 + IKH11 * Ppred11;

  P00 = newP00;
  P01 = newP01;
  P10 = newP10;
  P11 = newP11;
}

void printCSVHeader() {
  Serial.println("time_s,busVoltage_V,loadVoltage_V,current_A,power_W,soc_percent,vrc_V,ocv_V,critical_on,medium_on,deferrable_on");
}

void printCSV(float timeS, float busVoltageV, float loadVoltageV, float currentA, float powerW) {
  Serial.print(timeS, 3); Serial.print(",");
  Serial.print(busVoltageV, 4); Serial.print(",");
  Serial.print(loadVoltageV, 4); Serial.print(",");
  Serial.print(currentA, 5); Serial.print(",");
  Serial.print(powerW, 4); Serial.print(",");
  Serial.print(soc * 100.0, 2); Serial.print(",");
  Serial.print(vrc, 5); Serial.print(",");
  Serial.print(getOCV(soc), 4); Serial.print(",");
  Serial.print(digitalRead(MOSFET_CRITICAL)); Serial.print(",");
  Serial.print(digitalRead(MOSFET_MEDIUM)); Serial.print(",");
  Serial.println(digitalRead(MOSFET_DEFERRABLE));
}

void printSerialPlotter(float loadVoltageV, float currentA) {
  // Good for Arduino IDE Serial Plotter.
  Serial.print("SOC_percent:");
  Serial.print(soc * 100.0, 2);
  Serial.print("\tVoltage_V:");
  Serial.print(loadVoltageV, 4);
  Serial.print("\tCurrent_A:");
  Serial.print(currentA, 4);
  Serial.print("\tVRC_V:");
  Serial.println(vrc, 5);
}

void setup() {
  Serial.begin(115200);
  delay(1000);

  Wire.begin(SDA_PIN, SCL_PIN);

  pinMode(MOSFET_CRITICAL, OUTPUT);
  pinMode(MOSFET_MEDIUM, OUTPUT);
  pinMode(MOSFET_DEFERRABLE, OUTPUT);

  digitalWrite(MOSFET_CRITICAL, LOW);
  digitalWrite(MOSFET_MEDIUM, LOW);
  digitalWrite(MOSFET_DEFERRABLE, LOW);

  if (!ina226.init()) {
    Serial.println("INA226 not found. Check wiring, address, SDA, and SCL.");
    while (1);
  }

  ina226.setAverage(INA226_AVERAGE_16);
  ina226.setConversionTime(INA226_CONV_TIME_1100);
  ina226.setMeasureMode(INA226_CONTINUOUS);

  // Change 0.1 if your INA226 board uses a different shunt resistor.
  ina226.setResistorRange(0.1, 1.3);

  startTime = millis();
  lastEKFTime = millis();
  lastLogTime = millis();

  if (!SERIAL_PLOTTER_MODE) {
    Serial.println("INA226 connected. Starting EKF CSV logging.");
    printCSVHeader();
  } else {
    Serial.println("SOC_percent:0\tVoltage_V:0\tCurrent_A:0\tVRC_V:0");
  }
}

void loop() {
  unsigned long now = millis();

  if (now - lastEKFTime >= EKF_INTERVAL_MS) {
    float dt = (now - lastEKFTime) / 1000.0;
    lastEKFTime = now;

    // Load shedding based on current SOC estimate.
    updateLoads(soc);

    // -------------------- INA226 readings --------------------
    float shuntVoltage_mV = ina226.getShuntVoltage_mV();
    float busVoltage_V = ina226.getBusVoltage_V();
    float current_mA = ina226.getCurrent_mA();
    float power_mW = ina226.getBusPower();

    float current_A = CURRENT_SIGN * (current_mA / 1000.0);
    float power_W = power_mW / 1000.0;

    // Same voltage calculation from your starter code.
    float loadVoltage_V = busVoltage_V + (shuntVoltage_mV / 1000.0);

    // Run EKF using measured voltage and current.
    runEKF(loadVoltage_V, current_A, dt);


    printSerialPlotter(loadVoltage_V, current_A);
    
  /*else if (now - lastLogTime >= LOG_INTERVAL_MS) {
      lastLogTime = now;
      float timeS = (now - startTime) / 1000.0;
      printCSV(timeS, busVoltage_V, loadVoltage_V, current_A, power_W);
  }
      */
    
  }
}
