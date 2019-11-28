////
////  integration.cpp
////  nmr
////
////  Created by Eduardo Rubio on 10/17/19.
////  Copyright © 2019 Eduardo Rubio. All rights reserved.
////
//
#include <stdio.h>

/**
 *Integration Methods:
 * The following four functions are four distinct numerical methods for integration
 * The functions are predicated on using a pre-calculated cubic spline interpolant
 */
double adaptiveQuadrature( double x0, double x1, int splindex ){
    double app  = 0.0;
    double TOL = options.tol;
    int n = 50;
    int i = 0;
    double tol[n],a[n],h[n],fa[n],fc[n],fb[n],s[n],l[n];
    double v1,v2,v3,v4,v5,v6,v7,v8,s1,s2,fd,fe;
    // Step 1
    tol[i] = 10 * TOL;
    a[i] = x0;
    h[i] = (x1-x0)/2;
    fa[i] = f( x0 - _data.x[ splindex ],
              _data.y[ splindex ] - options.baseline,
              splineCo.B[ splindex ],
              splineCo.C[ splindex ],
              splineCo.D[ splindex ]
              );
    fc[i] = f(x0 + h[i] - _data.x[ splindex ],
              _data.y[ splindex ] - options.baseline,
              splineCo.B[ splindex ],
              splineCo.C[ splindex ],
              splineCo.D[ splindex ]
              );
    fb[i] = f( x1 - _data.x[ splindex ],
              _data.y[ splindex ] - options.baseline,
              splineCo.B[ splindex ],
              splineCo.C[ splindex ],
              splineCo.D[ splindex ]);
    s[i] = h[i]*(fa[i] + 4*fc[i] + fb[i])/3;
    l[i] = 1;
    // Step 2
    while (i>-1){
        // Step 3
        fd = f((a[i] + h[i]/2)- _data.x[ splindex],
               _data.y[ splindex ] - options.baseline,
               splineCo.B[ splindex ],
               splineCo.C[ splindex ],
               splineCo.D[ splindex ]);
        fe = f((a[i] + 3*h[i]/2) - _data.x[ splindex ],
               _data.y[ splindex ] - options.baseline,
               splineCo.B[ splindex ],
               splineCo.C[ splindex ],
               splineCo.D[ splindex ]);
        s1 = h[i]*(fa[i] + 4*fd + fc[i])/6;
        s2 = h[i]*(fc[i] + 4*fe + fb[i])/6;
        v1 = a[i];
        v2 = fa[i];
        v3 = fc[i];
        v4 = fb[i];
        v5 = h[i];
        v6 = tol[i];
        v7 = s[i];
        v8 = l[i];
        // Step 4
        i--;
        //Step 5
        if (abs(s1+s2-v7)<v6)
            app += s1+s2;
        else
            if (v8>=n)
                std::cout << "Level Exceeded" << std::endl;
            else{
                i++;
                a[i] = v1+v5;
                fa[i] = v3;
                fc[i] = fe;
                fb[i] = v4;
                h[i] = v5/2;
                tol[i] = v6/2;
                s[i] = s2;
                l[i] = v8+1;
                
                i++;
                a[i] = v1;
                fa[i] = v2;
                fc[i] = fd;
                fb[i] = v3;
                h[i] = h[i-1];
                tol[i] = tol[i-1];
                s[i] = s1;
                l[i] = l[i-1];
            }
    }
    //Step 6
    return app;
}

double romberg( double a, double b, int splindex ){
    int n = 10;
    double tol = options.tol;
    double h  = b - a;
    std::vector<std::vector<double> > R( n , std::vector<double> (n, 0));
    // Step 1
    double fa = f( a - _data.x[ splindex ],
                  _data.y[ splindex ] - options.baseline,
                  splineCo.B[ splindex ],
                  splineCo.C[ splindex ],
                  splineCo.D[ splindex ]
                  );
    double fb = f( b - _data.x[ splindex ],
                  _data.y[ splindex ]- options.baseline,
                  splineCo.B[ splindex ],
                  splineCo.C[ splindex ],
                  splineCo.D[ splindex ]
                  );
    R[0][0] = (h/2) * (fa+fb);
    
    // Step 3
    for (int i = 1; i<n;i++) {
        // Step 4
        double sum = 0;
        for( int k = 1; k <= (pow(2, i-1)); k++ )
            sum += f( (a + (k-0.5)*h) - _data.x[ splindex ],
                     _data.y[ splindex ]- options.baseline,
                     splineCo.B[ splindex ],
                     splineCo.C[ splindex ],
                     splineCo.D[ splindex ]
                     );
        R[1][0] = 0.5 * ( R[0][0] + h * sum );
        // Step 5
        for( int j = 1; j <= i; j++ )
            R[1][j] = R[1][j-1] + (R[1][j-1]-R[0][j-1])/(pow(4.0,j)-1);
        // Step 7
        h = h / 2;
        //Step 8
        for (int j = 0; j <=i; j++)
            R[0][j] = R[1][j];
        if( fabs(R[i][i] - R[i-1][i-1]) < tol ) return R[1][i];
    }
    return -1;
}

double compositeSimsons( double a, double b, int splindex ){
    int n = 10;
    double h = (b-a)/n;
    double x, xi;
    double xi0 = f(a - _data.x[ splindex ],
                   _data.y[ splindex ] - options.baseline,
                   splineCo.B[ splindex ],
                   splineCo.C[ splindex ],
                   splineCo.D[ splindex ]
                   )
    + f(b - _data.x[ splindex ],
        _data.y[ splindex ] - options.baseline,
        splineCo.B[ splindex ],
        splineCo.C[ splindex ],
        splineCo.D[ splindex ]
        );
    double xi1 = 0;
    double xi2 = 0;
    for( int i = 1; i <= (n-1); i++ ){
        x = a + i*h;
        if( i % 2 == 0 ) xi2 = xi2 + f( x - _data.x[ splindex ],
                                       _data.y[ splindex ] - options.baseline,
                                       splineCo.B[ splindex ],
                                       splineCo.C[ splindex ],
                                       splineCo.D[ splindex ]
                                       );
        else xi1 = xi1 + f(x - _data.x[ splindex ],
                           _data.y[ splindex ] - options.baseline,
                           splineCo.B[ splindex ],
                           splineCo.C[ splindex ],
                           splineCo.D[ splindex ]
                           );
    }
    xi = h*(xi0 + 2*xi2 + 4*xi1);
    xi = xi/3;
    return xi;
}

double gaussLaguerre( double a, double b, int splindex ){
    double X15[15] = {
        9.33078120172818E-02,
        4.92691740301884E-01,
        1.21559541207095E+00,
        2.26994952620374E+00,
        3.66762272175144E+00,
        5.42533662741355E+00,
        7.56591622661307E+00,
        1.01202285680191E+01,
        1.31302824821757E+01,
        1.66544077083300E+01,
        2.07764788994488E+01,
        2.56238942267288E+01,
        3.14075191697539E+01,
        3.85306833064860E+01,
        4.80260855726858E+01
    };
    double W15[15] = {
        2.18234885940085E-01,
        3.42210177922883E-01,
        2.63027577941712E-01,
        1.26425818105934E-01,
        4.02068649210009E-02,
        8.56387780361183E-03,
        1.21243614721425E-03,
        1.11674392344252E-04,
        6.45992676202291E-06,
        2.22631690709627E-07,
        4.22743038497937E-09,
        3.92189726704109E-11,
        1.45651526407312E-13,
        1.48302705111330E-16,
        1.60059490621113E-20
    };
    double sum = 0.0;
    for( int i = 0; i <= 15; i++ ){
        sum += W15[i] * f(((b-a)*X15[i]+(b+a))/2.0 - _data.x[ splindex ],
                          _data.y[ splindex ] - options.baseline,
                          splineCo.B[ splindex ],
                          splineCo.C[ splindex ],
                          splineCo.D[ splindex ]
                          ) * (b-a)/2.0;
    }
    return sum;
}

/**
 * Iterate through all the peaks in the Peak vector and calculate
 * their peak sums through the user designated method of intergration
 *
 * The method of integration is user-defined, but each method will iterate through
 * the encompassing cubic spline interpolants the peak traverses
 *
 * Integration methods include:
 *      1 ** Adaptive Quadrature
 *      2 ** Romberg's Method
 *      3 ** Gaussian Laguerre Integration
 *      4 ** Composite Simpson's Method
 */
void integratePeaks( ){
    for( int peakID = 0; peakID < peaks.size(); peakID++ ){
        double sum = 0;
        Peak tmp = peaks[ peakID ];
        // If peak bound by domain
        if( !tmp.isComplete ){
            tmp.indexB = _data.x.size();
            tmp.rootB = _data.x[ tmp.indexB ];
        }
        switch( options.integrationType ) {
                
            case 0: //Adaptive Quadrature
                sum += adaptiveQuadrature(  tmp.rootA, _data.x[ tmp.indexA + 1 ], tmp.indexA );
                for( int j = tmp.indexA + 1; j < tmp.indexB; j++ )
                    sum += adaptiveQuadrature(_data.x[ j ], _data.x[ j + 1 ], j );
                sum += adaptiveQuadrature( _data.x[ tmp.indexB ], tmp.rootB, tmp.indexB );
                break;
                
            case 1: // Romberg
                sum += romberg(  tmp.rootA, _data.x[ tmp.indexA + 1 ], tmp.indexA );
                for( int j = tmp.indexA + 1; j < tmp.indexB; j++ )
                    sum += romberg(_data.x[ j ], _data.x[ j + 1 ], j );
                sum += romberg( _data.x[ tmp.indexB ], tmp.rootB, tmp.indexB );
                break;
                
            case 2: //Gaussian Laguerre
                sum += gaussLaguerre(  tmp.rootA, _data.x[ tmp.indexA + 1 ], tmp.indexA );
                for( int j = tmp.indexA + 1; j < tmp.indexB; j++ )
                    sum += gaussLaguerre(_data.x[ j ], _data.x[ j + 1 ], j );
                sum += gaussLaguerre( _data.x[ tmp.indexB ], tmp.rootB, tmp.indexB );
                break;
                
            case 3: //Composite Simpson's
                sum += compositeSimsons(  tmp.rootA, _data.x[ tmp.indexA + 1 ], tmp.indexA );
                for( int j = tmp.indexA + 1; j < tmp.indexB; j++ )
                    sum += compositeSimsons(_data.x[ j ], _data.x[ j + 1 ], j );
                sum += compositeSimsons( _data.x[ tmp.indexB ], tmp.rootB, tmp.indexB );
                break;
                
            default:
                printf("Error: Undefined integration options requested\n");
                break;
        }
        peaks[ peakID ].manifold = fabs( sum );
    }
}

/**
 * Iterate through the now fully calculated peaks and find hydrogen values
 */
void findHydrogens( ){
    double min_area = __DBL_MAX__;
    for( int i = 0; i < peaks.size(); i++ )
        min_area = fmin( min_area, peaks[i].manifold);
    for( int i = 0; i < peaks.size(); i++ )
        peaks[i].hydrogens = peaks[i].manifold / min_area;
}
















