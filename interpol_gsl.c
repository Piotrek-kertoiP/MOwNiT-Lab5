#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <sys/resource.h>
#include <unistd.h>

#define NodesAmount 50
#define InterNodeGap 5
#define REPETITIONS 10
#define VERBOSE 1
long double lagrangeCoefficients [NodesAmount];
double newtonCoefficients  [NodesAmount ];
FILE * resultFile;
FILE * timesFile;
//-------------------------------------------------------------------------------------NODES GENERATOR------------------------------------------------------------------------------------
double doubleRandRange(double minN, double maxN){
    return (double)rand()/RAND_MAX * (maxN - minN) + minN;
}

void randNodes(double* x ,double* y ){
    srand(time(NULL));
    x[0] = 0.0;
    y[0] = doubleRandRange(-100, 100);
    for ( int i = 0; i < NodesAmount - 1; i++)
    {
        x[i+1] = x[i] + InterNodeGap;
        y[i+1] = doubleRandRange(-100, 100);
    }

}
//-------------------------------------------------------------------------------------GSL INTERPOLATION----------------------------------------------------------------------------------
void gslInterpolation(double* x, double* y, int stepsPerGap){
    double tmpX, tmpY, h;
    // interpolation
    {
        gsl_interp_accel *acc = gsl_interp_accel_alloc ();
        gsl_interp * workspace = gsl_interp_alloc(gsl_interp_polynomial, NodesAmount);

        gsl_interp_init(workspace, x, y, NodesAmount);
        h = (x[1]-x[0]) / stepsPerGap;
        for(int i = 0; i <= (stepsPerGap *(NodesAmount-1)); i++)
        {
            tmpX = x[0] + h * i;
            tmpY = gsl_interp_eval(workspace, x, y, tmpX, acc);
            if( VERBOSE )
                printf("GSL,%lf,%lf \n", tmpX, tmpY);
            fprintf(resultFile,"GSL,%lf,%lf \n", tmpX, tmpY);
        }

        gsl_interp_free (workspace);
        gsl_interp_accel_free (acc);
    }
}
//------------------------------------------------------------------------------------SPLINE INTERPOLATION-------------------------------------------------------------------------------
void gslSplineInterpolation(double* x, double* y, int stepsPerGap ){
    double tmpX, tmpY, h;
    // interpolation
    {
        gsl_interp_accel *acc = gsl_interp_accel_alloc ();
        gsl_spline * workspace = gsl_spline_alloc(gsl_interp_cspline, NodesAmount);

        gsl_spline_init(workspace, x, y, NodesAmount);
        h = (x[1]-x[0])/stepsPerGap;
        for(int i = 0; i <= (stepsPerGap*(NodesAmount-1)); i++)
        {
            tmpX = x[0] + h * i;
            tmpY = gsl_spline_eval(workspace, tmpX, acc);
            if( VERBOSE )
                printf("Spline,%lf,%lf \n", tmpX, tmpY);
            fprintf(resultFile,"Spline,%lf,%lf \n", tmpX, tmpY);
        }

        gsl_spline_free (workspace);
        gsl_interp_accel_free (acc);
    }
}
//-----------------------------------------------------------------------------------LAGRANGE INTERPOLATION-------------------------------------------------------------------------------
void getLagrangeCoefficients(double* x, double* y){
    // c[] are the coefficients of P(x) = (x-x[0])(x-x[1])...(x-x[n-1])
    long double c[NodesAmount+1];
    for (int i = 0; i < NodesAmount; i++){
        lagrangeCoefficients[i] = 0.0;
    }
    c[0] = 1.0;
    for ( int i = 0; i < NodesAmount; i++){
        for ( int j = i; j > 0; j--){
            c[j] = c[j-1] - c[j] * x[i];
        }
        c[0] *= -x[i];
        c[i+1] = 1;
    }
    long double tc[NodesAmount+1];
    for ( int i = 0; i < NodesAmount; i++) {
        // d = (x[i]-x[0])...(x[i]-x[i-1])(x[i]-x[i+1])...(x[i]-x[n-1])
        long double d = 1.0;
        for ( int j = 0; j < NodesAmount; j++) {
            if (i != j) {
                d *= (x[i] - x[j]);
            }
        }
        if (d == 0.0) {
            // This happens only when two abscissas are identical.
            for ( int k = 0; k < NodesAmount; k++) {
                if ((i != k) && (x[i] == x[k])) {
                    printf("LocalizedFormats.IDENTICAL_ABSCISSAS_DIVISION_BY_ZERO,%d %d %lf",i, k, x[i]);
                    exit(-1);
                }

            }
        }

        long double t = y[i] / d;
        // Lagrange polynomial is the sum of n terms, each of which is a
        // polynomial of degree n-1. tc[] are the coefficients of the i-th
        // numerator Pi(x) = (x-x[0])...(x-x[i-1])(x-x[i+1])...(x-x[n-1]).
        tc[NodesAmount-1] = c[NodesAmount];     // actually c[n] = 1
        lagrangeCoefficients[NodesAmount-1] += t * tc[NodesAmount-1];
        for ( int j = NodesAmount-2; j >= 0; j--) {
            tc[j] = c[j+1] + tc[j+1] * x[i];
            lagrangeCoefficients[j] += (t * tc[j]);
        }
    }
    return;
}
double lagrangePolynomialValue(long double* coefficients, int numberOfCoefficients, long double x)
{
    long double result = 0.0;

    for(int i = numberOfCoefficients - 1; i >= 0; i--)
    {
        result = result * x + coefficients[i];
    }
    return result;
}

void lagrangeInterpolation(double x[],double y[],int stepsPerGap)
{
    long double tmpX, tmpY, h;
    h = (x[1]-x[0]) / stepsPerGap;
    for(int i = 0; i <= ( stepsPerGap * (NodesAmount-1) ); i++)
    {
        tmpX = x[0] + h * i;
        tmpY = lagrangePolynomialValue(lagrangeCoefficients, NodesAmount, tmpX);
        if( VERBOSE )
            printf("Lagrange,%Lf,%Lf \n", tmpX, tmpY);
        fprintf(resultFile,"Lagrange,%Lf,%Lf \n", tmpX, tmpY);
    }
}
//-----------------------------------------------------------------------------------NEWTON INTERPOLATION--------------------------------------------------------------------------------
void getNewtonCoefficients(double* x, double* y){
    double newtonCoeffs [NodesAmount-1][NodesAmount-1];

    // zero table
    for (int i = 0; i < (NodesAmount-1); i++){
        for (int j = 0; j < (NodesAmount-1); j++){
            newtonCoeffs[i][j] = 0;
        }
    }

    // first columns
    for (int i = 0; i < (NodesAmount-1); i++){
        newtonCoeffs[0][i] = (y[i+1] - y[i]) / (x[i+1] - x[i]);
    }

    for (int i = 1; i < (NodesAmount-1); i++){
        int counter = 0;
        for(int j = (NodesAmount-2-i); j >= 0; j--){
            newtonCoeffs[i][j] = (newtonCoeffs[i-1][j+1] - newtonCoeffs[i-1][j]) / (x[(NodesAmount-1-counter)] - x[(NodesAmount-2-i-counter)]);
            counter++;
        }
    }

    newtonCoefficients[0] = y[0];
    for (int i = 1; i < NodesAmount; i++){
        newtonCoefficients[i] = newtonCoeffs[i-1][0];
    }
}
double newtonPolynomialValue(double* coefficients, double* nodesX, int summandsAmount, double x)
{
    double p = 1;
    double result = 0.0;
    for (int i = 0; i < summandsAmount; i++)
    {
        result += (p * coefficients[i]);
        p = p * (x - nodesX[i]);
    }
    return result;

}
void newtonInterpolation(double* x, double* y, int stepsPerGap)
{
    double tmpX, tmpY, h;
    h = (x[1]-x[0]) / stepsPerGap;
    for(int i = 0; i <= (stepsPerGap * (NodesAmount-1) ); i++)
    {
        tmpX = x[0] + h * i;
        tmpY = newtonPolynomialValue(newtonCoefficients, x, NodesAmount, tmpX);
        if( VERBOSE )
            printf("Newton,%lf,%lf \n", tmpX,tmpY);
        fprintf(resultFile,"Newton,%lf, %lf \n", tmpX,tmpY);
    }
}
//------------------------------------------------------------------------------------AKIMA INTERPOLATION--------------------------------------------------------------------------------
void gslAkimaInterpolation(double* x, double* y, int stepsPerGap){
    double tmpX, tmpY, h;
    // interpolation
    {
        gsl_interp_accel *acc = gsl_interp_accel_alloc ();
        gsl_interp * workspace = gsl_interp_alloc(gsl_interp_akima, NodesAmount);

        gsl_interp_init(workspace, x, y, NodesAmount);
        h = (x[1]-x[0]) / stepsPerGap;
        for(int i = 0; i <= (stepsPerGap*(NodesAmount-1)); i++){
            tmpX = x[0] + h * i;
            tmpY = gsl_interp_eval(workspace, x, y, tmpX, acc);
            if( VERBOSE )
                printf("Akima,%lf,%lf \n", tmpX,tmpY);
            fprintf(resultFile,"Akima,%lf,%lf \n", tmpX,tmpY);
        }

        gsl_interp_free (workspace);
        gsl_interp_accel_free (acc);
    }
}
//--------------------------------------------------------------------------------------TIME MEASUREMENT---------------------------------------------------------------------------------
struct timeval lagrangeTimeStart;
struct timeval lagrangeTimeEnd;

struct timeval newtonTimeStart;
struct timeval newtonTimeEnd;

struct timeval gslTimeStart;
struct timeval gslTimeEnd;

struct timeval splineTimeStart;
struct timeval splineTimeEnd;

struct timeval akimaTimeStart;
struct timeval akimaTimeEnd;

double convertToDouble( struct timeval start, struct timeval end){
    double elapsedTime = (end.tv_sec - start.tv_sec) * 1000.0;      // sec to ms
    elapsedTime += (end.tv_usec - start.tv_usec) / 1000.0;   // us to ms
    return elapsedTime;
}
//------------------------------------------------------------------------------------SAVING TIMES TO FILE-------------------------------------------------------------------------------
void saveTimesToFile(double lagrange, double newton, double gsl, double spline, double akima, int N){
    fprintf(timesFile, "%d, Lagrange, %f\n", N, lagrange);
    fprintf(timesFile, "%d, Newton, %f\n", N, newton);
    fprintf(timesFile, "%d, GSL, %f\n", N, gsl);
    fprintf(timesFile, "%d, Spline, %f\n", N, spline);
    fprintf(timesFile, "%d, Akima, %f\n", N, akima);
}
//-------------------------------------------------------------------------------------------MAIN----------------------------------------------------------------------------------------

int main (){
    if  ((resultFile = fopen ("result.txt", "w")) == NULL ){
        printf("Result file open error \n");
        exit(-1);
    }
    if( (timesFile = fopen ("times.txt", "a+")) == NULL){
        printf("Times file open error \n");
        exit(-1);
    }

    double x[NodesAmount];
    double y[NodesAmount];
    randNodes(x,y);

    if( VERBOSE )
        printf("alg,x,y \n");
    //fprintf(resultFile,"alg,x,y \n");
    for ( int i = 0; i < NodesAmount; i++)
    {
        if( VERBOSE )
            printf("Generated_nodes, %lf,%lf \n", x[i],y[i]);
        fprintf(resultFile,"Generated_nodes, %lf,%lf \n", x[i],y[i]);
    }
                
    for( int j = 0; j < REPETITIONS; j++){
        gettimeofday(&lagrangeTimeStart, NULL);
        getLagrangeCoefficients(x,y);
        lagrangeInterpolation(x,y,InterNodeGap);
        gettimeofday(&lagrangeTimeEnd, NULL);
        double lagrange = convertToDouble(lagrangeTimeStart, lagrangeTimeEnd);

        gettimeofday(&newtonTimeStart, NULL);
        getNewtonCoefficients(x,y);
        newtonInterpolation(x,y,InterNodeGap);
        gettimeofday(&newtonTimeEnd, NULL);
        double newton = convertToDouble(newtonTimeStart, newtonTimeEnd);

        gettimeofday(&gslTimeStart, NULL);
        gslInterpolation(x,y,InterNodeGap);
        gettimeofday(&gslTimeEnd, NULL);
        double gsl = convertToDouble(gslTimeStart, gslTimeEnd);

        gettimeofday(&splineTimeStart, NULL);
        gslSplineInterpolation(x,y,InterNodeGap);
        gettimeofday(&splineTimeEnd, NULL);
        double spline = convertToDouble(splineTimeStart, splineTimeEnd);

        gettimeofday(&akimaTimeStart, NULL);
        gslAkimaInterpolation(x,y,InterNodeGap);
        gettimeofday(&akimaTimeEnd, NULL);
        double akima = convertToDouble(akimaTimeStart, akimaTimeEnd);

        saveTimesToFile(lagrange, newton, gsl, spline, akima, NodesAmount);
    }
    return 0;
}
