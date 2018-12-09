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
double gsl_coefficients [NodesAmount];
FILE * resultFile;
FILE * times_file;

double double_rand_range(double min_n, double max_n)
{
    return (double)rand()/RAND_MAX * (max_n - min_n) + min_n;
}

void rand_nodes(double x [] ,double y [] )
{
    int i;
    srand(time(NULL));
    x[0]=0;
    y[0]=double_rand_range(-100,100);
    for ( i = 0; i < NodesAmount - 1; i++)
    {
        x[i+1] = x[i] + InterNodeGap;
        y[i+1] = double_rand_range(-100,100);
    }

}

void gsl_polym_interpol(double x [] ,double y [],int step )
{
    int i;
    double h;
    double xtmp,ytmp;
    // interpolation
    {
        gsl_interp_accel *acc = gsl_interp_accel_alloc ();
        gsl_interp * workspace = gsl_interp_alloc(gsl_interp_polynomial, NodesAmount);

        gsl_interp_init(workspace, x, y, NodesAmount);
        h = (x[1]-x[0]) / step;
        for(i = 0; i <= (step *(NodesAmount-1)); i++)
        {
            xtmp = x[0] + h *i;
            ytmp = gsl_interp_eval(workspace, x, y, xtmp, acc);
            if( VERBOSE )
                printf("GSL,%lf,%lf \n", xtmp,ytmp);
            fprintf(resultFile,"GSL,%lf,%lf \n", xtmp,ytmp);
        }

        gsl_interp_free (workspace);
        gsl_interp_accel_free (acc);
    }
}

void gsl_spline_interpol(double x [] ,double y [],int step )
{
    int i;
    double h;
    double xtmp,ytmp;
    // interpolation
    {
        gsl_interp_accel *acc = gsl_interp_accel_alloc ();
        gsl_spline * workspace = gsl_spline_alloc(gsl_interp_cspline, NodesAmount);

        gsl_spline_init(workspace, x, y, NodesAmount);
        h = (x[1]-x[0])/step;
        for(i = 0; i <= (step*(NodesAmount-1)); i++)
        {
            xtmp = x[0] + h *i;
            ytmp = gsl_spline_eval(workspace, xtmp, acc);
            if( VERBOSE )
                printf("Spline,%lf,%lf \n", xtmp,ytmp);
            fprintf(resultFile,"Spline,%lf,%lf \n", xtmp,ytmp);
        }

        gsl_spline_free (workspace);
        gsl_interp_accel_free (acc);
    }
}

// LaGrange init
void getLagrangeCoeffs(double x[],double y[])  {
    int i,j,k;
    // c[] are the coefficients of P(x) = (x-x[0])(x-x[1])...(x-x[n-1])
    long double c[NodesAmount+1];
    for (i=0;i<NodesAmount;i++)
    {
        lagrangeCoefficients[i] = 0.0;
    }
    c[0] = 1.0;
    for ( i = 0; i < NodesAmount; i++)
    {
        for ( j = i; j > 0; j--)
        {
            c[j] = c[j-1] - c[j] * x[i];
        }
        c[0] *= -x[i];
        c[i+1] = 1;
    }
    long double tc[NodesAmount+1];
    for ( i = 0; i < NodesAmount; i++) {
        // d = (x[i]-x[0])...(x[i]-x[i-1])(x[i]-x[i+1])...(x[i]-x[n-1])
        long double d = 1.0;
        for ( j = 0; j < NodesAmount; j++) {
            if (i != j) {
                d *= (x[i] - x[j]);
            }
        }
        if (d == 0.0) {
            // This happens only when two abscissas are identical.
            for ( k = 0; k < NodesAmount; ++k) {
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
        for ( j = NodesAmount-2; j >= 0; j--) {
            tc[j] = c[j+1] + tc[j+1] * x[i];
            lagrangeCoefficients[j] += (t * tc[j]);
        }
    }
    return;
}

void getNetwonCoeffs(double x[],double y[])  {
    //  ilorazy roznicowe nie cala jest potrzebna, tylko  polowa
    double newton_coeffs [NodesAmount-1][NodesAmount-1];
    int i,j;

    // zero table
    for (i=0;i<(NodesAmount-1);i++)
        for (j=0;j<(NodesAmount-1);j++)
            newton_coeffs[i][j]=0;

    // first columns
    for (i=0;i<(NodesAmount-1);i++)
        newton_coeffs[0][i]=(y[i+1]-y[i])/(x[i+1]-x[i]);


    for (i=1;i<(NodesAmount-1);i++)
    {
        int counter=0;
        for(j=(NodesAmount-2-i);j>=0;j--)
        {
            newton_coeffs[i][j]=(newton_coeffs[i-1][j+1]-newton_coeffs[i-1][j])/(x[(NodesAmount-1-counter)]-x[(NodesAmount-2-i-counter)]  );
            counter++;
        }
    }

    newtonCoefficients[ 0  ]=y[0];
    for (i=1;i<NodesAmount;i++)
    {
        newtonCoefficients[ i ]=newton_coeffs[i-1][0];
    }
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

void lagrangeInterpol(double x[],double y[],int stepsPerGap)
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
void newtonInterpol(double* x, double* y, int stepsPerGap)
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

void gsl_akima_interpol(double x [] ,double y [],int step )
{
    int i;
    double h;
    double xtmp,ytmp;
    // interpolation
    {
        gsl_interp_accel *acc = gsl_interp_accel_alloc ();
        gsl_interp * workspace = gsl_interp_alloc(gsl_interp_akima, NodesAmount);

        gsl_interp_init(workspace, x, y, NodesAmount);
        h = (x[1]-x[0])/step;
        for(i = 0; i <= (step*(NodesAmount-1)); i++)
        {
            xtmp = x[0] + h *i;
            ytmp = gsl_interp_eval(workspace, x, y, xtmp, acc);
            if( VERBOSE )
                printf("Akima,%lf,%lf \n", xtmp,ytmp);
            fprintf(resultFile,"Akima,%lf,%lf \n", xtmp,ytmp);
        }

        gsl_interp_free (workspace);
        gsl_interp_accel_free (acc);
    }

}


//------------------------------------------TIME MEASUREMENT--------------------
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
//-------------------------------------------------------------------------SAVING TIMES TO FILE
void saveTimesToFile(double lagrange, double newton, double gsl, double spline, double akima, int N){
    fprintf(times_file, "%d, Lagrange, %f\n", N, lagrange);
    fprintf(times_file, "%d, Newton, %f\n", N, newton);
    fprintf(times_file, "%d, GSL, %f\n", N, gsl);
    fprintf(times_file, "%d, Spline, %f\n", N, spline);
    fprintf(times_file, "%d, Akima, %f\n", N, akima);
}
//-----------------------------------------------------------------------------

int main ()
{
    if  ((resultFile = fopen ("result.txt", "w")) == NULL )
    {
        printf("Result file open error \n");
        exit(-1);
    }
    if( (times_file = fopen ("times.txt", "a+")) == NULL)
    {
        printf("Times file open error \n");
        exit(-1);
    }

    double x[NodesAmount];
    double y[NodesAmount];
    rand_nodes(x,y);

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
        getLagrangeCoeffs(x,y);
        lagrangeInterpol(x,y,InterNodeGap);
        gettimeofday(&lagrangeTimeEnd, NULL);
        double lagrange = convertToDouble(lagrangeTimeStart, lagrangeTimeEnd);

        gettimeofday(&newtonTimeStart, NULL);
        getNetwonCoeffs(x,y);
        newtonInterpol(x,y,InterNodeGap);
        gettimeofday(&newtonTimeEnd, NULL);
        double newton = convertToDouble(newtonTimeStart, newtonTimeEnd);

        gettimeofday(&gslTimeStart, NULL);
        gsl_polym_interpol(x,y,InterNodeGap);
        gettimeofday(&gslTimeEnd, NULL);
        double gsl = convertToDouble(gslTimeStart, gslTimeEnd);

        gettimeofday(&splineTimeStart, NULL);
        gsl_spline_interpol(x,y,InterNodeGap);
        gettimeofday(&splineTimeEnd, NULL);
        double spline = convertToDouble(splineTimeStart, splineTimeEnd);

        gettimeofday(&akimaTimeStart, NULL);
        gsl_akima_interpol(x,y,InterNodeGap);
        gettimeofday(&akimaTimeEnd, NULL);
        double akima = convertToDouble(akimaTimeStart, akimaTimeEnd);

        saveTimesToFile(lagrange, newton, gsl, spline, akima, NodesAmount);
    }
    return 0;
}
