#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// to generate more cycle, increase cycle number
// then duplicate the input based on cycle number
// for example : in_a[] = {Va[0], .... , Va[DATA_LENGTH - 1],
// Va[0], .... , Va[DATA_LENGTH - 1],
// . . .};

#define CYCLE 1
#define DATA_LENGTH 20
#define PI 3.14159265359

float Va[] = {
    156.63, 246.59, 294.72, 305.51, 300.66,
    268.03, 204.18, 125.41, 42.954, -48.322,
    -154.08, -243.95, -293.12, -303.09, -297.98,
    -264.13, -202.1, -122.25, -39.893, 51.818
};

float Vb[] = {
    -308.4, -280.19, -240.66, -186.6, -99.744,
    -0.54547, 92.853, 181.46, 262.05, 312.39,
    311.44, 283.76, 245.04, 188.62, 102.16,
    2.9662, -89.395, -176.17, -259.16, -309.96
};

float Vc[] = {
    156.11, 82.694, -21.783, -128.37, -213.06,
    -269.49, -309.58, -313.4, -273.73, -214.81,
    -154.29, -79.64, 24.679, 132.16, 216.63,
    274.14, 311.11, 315.76, 276.27, 216.22
};

typedef struct _DDATA{
    float *in_a;
    float *in_b;
    float *in_c;
    float *F_est;
    float *Theta_est;
    float *Harmonics;
    float Ts;
    float Kc1; // Kc are controller gains
    float Kc2; // choose your controller and
    float Kc3; // gains accordingly to get satisfied result
}DDATA;

DDATA ddata = {
    // Removed incorrect pointer dererence since array name itself points to its first element
    .in_a = Va,
    .in_b = Vb,
    .in_c = Vc,
    .Ts = 0.001,
};

// Function estimateFrequencyAndTheta calculates estimates
// for Theta_est and F_est on each data instance.
// Initially, the estimates may not be very accurate,
// but they shall improve with more input data.
// Thus, the size of the Estimated Theta and Freq arrays
// should match the input data size (DATALENGTH * CYCLE)
void estimateFrequencyAndTheta(DDATA *d, int dataSize){

    // Use dynamic allocation for F_est & Theta_est then depending on actual target
    // consider posibility of using static allocation
    d->F_est = (float *)malloc(dataSize * sizeof(float));
    d->Theta_est = (float *)malloc(dataSize * sizeof(float));

    float alpha, beta;

    for (int i = 0; i < dataSize; i++) {

        // Using Clarke transformation to get alpha and beta using values a, b & c values
        alpha = (2.0/3.0) * (d->in_a[i] - 0.5 * (d->in_b[i] + d->in_c[i]));
        beta = (2.0/3.0) * ((sqrt(3)/2.0) * (d->in_b[i] - d->in_c[i]));

        // Calculate phase angle using atan2(beta, alpha)
        d->Theta_est[i] = atan2(beta, alpha);

        // Unwrap phase angles to ensure continuity
        if (i == 0){
          // First sample keep as it is, no Unwrap needed
        } else if (d->Theta_est[i] - d->Theta_est[i - 1] > PI) {
            d->Theta_est[i] -= 2 * PI;
        } else if (d->Theta_est[i] - d->Theta_est[i - 1] < -PI) {
            d->Theta_est[i] += 2 * PI;
        }

        if (i > 0){
            // Estimate frequency as dTheta/dt
            d->F_est[i] = (d->Theta_est[i] - d->Theta_est[i - 1]) / (2.0 * PI * d->Ts);
        } else{
            d->F_est[i] = 0;
        }
    }
}

// Function getHarmonicAmplitudes calculates the amplitudes of the 1st to 5th
// The output should consist of 5 data points, each representing the amplitud
// Additionally, you can use these 5 data points to calculate the Total Harmo
void getHarmonicAmplitudes(DDATA *d, int dataSize){
    // Implementation for getting harmonic amplitudes

}
int main()
{
    int i = 0;

    // for(i = 0; i < DATA_LENGTH * CYCLE; i++)
    {
        // Removed outside iteration since estimateFrequencyAndTheta itself iterates all sample points 
        estimateFrequencyAndTheta(&ddata, DATA_LENGTH * CYCLE);
    }

    // Print the Estimated Frequency by getting the last estimate value
    printf("Estimated Frequency (Hz): %.2f\n\n", ddata.F_est[(DATA_LENGTH * CYCLE)-1]);

    // Print the Estimated Phase angle over time
    printf("Estimated Phase angle (radians):\n");
    for (int i = 0; i < DATA_LENGTH * CYCLE; i++) {
        printf("%.2f, ", ddata.Theta_est[i]);
    }

    getHarmonicAmplitudes(&ddata, DATA_LENGTH * CYCLE);

    // Free allocated memory for F_est & Theta_est
    free(ddata.F_est);
    free(ddata.Theta_est);

    return 0;
}