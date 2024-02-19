#include "c74_min.h"
#include "../shared/signal_routing_objects.h"
#include <math.h>
#include <armadillo>
arma::fvec test(4, arma::fill::zeros);
class ladder : public object<ladder>, public sample_operator<4, 1>
{
private:
    arma::fvec::fixed<5> y;
    arma::fvec::fixed<5> y_est;
    arma::fvec::fixed<5> delta_y;
    arma::fvec::fixed<5> F;
    arma::fmat::fixed<5,5> J;
    
public:
    MIN_DESCRIPTION {"Transistor ladder filter with feedback."};
    MIN_TAGS {"filter, ladder, lpf"};
    MIN_AUTHOR {"Peter Corbett (Osprey Instruments)"};
    MIN_RELATED {"lores~, svf~, onepole~"};


    inlet<>  x {this, "(signal) Signal Input"};
    inlet<>  cutoff {this, "(signal) Cutoff Frequency"};
    inlet<>  r {this, "(signal) Resonance"};
    inlet<>  fbk {this, "(number) Feedback"};
    outlet<> output {this, "(signal) Output", "signal"};
        
    attribute<float> b { this, "bias", 0.65, range {-1.0, 1.0} };

    /*
    message<> dspsetup { this, "dspsetup",
        MIN_FUNCTION {
            number samplerate = args[0];
            //int vectorsize = args[1];
            //T = 1.0 / samplerate;
            return {};
        }
    };
    */
    
    // highpass filter
    float Fc_hp = 8.f;
    float g_hp = Fc_hp * M_PI / samplerate();
    
    // initial values
    /*
    fvec y(5);
    fvec y_est(5);
    fvec delta_y(5);
    fvec F(5);
    fmat J(5, 5);
     */
    float errorThresh = 0.000001; // max NR error (stopping condition)
        
    // accumulator states (capacitors)
    float s1 = 0.f;
    float s2 = 0.f;
    float s3 = 0.f;
    float s4 = 0.f;
    float s5 = 0.f;
        
        
    // Call operator: process a single sample
    sample operator()(sample x, sample cutoff, sample r, sample fbk)
    {
        float g = cutoff * M_PI / samplerate();
        
        int i = 0;
        int residue = 100;
        while(abs(residue) > errorThresh && i < 50)
        {
            // pre-compute tanh functions
            float tanh_x = tanh(x - r*y(3) + y(4));
            float tanh_y1 = tanh(y(0));
            float tanh_y2 = tanh(y(1));
            float tanh_y3 = tanh(y(2));
            float tanh_y4 = tanh(y(3));
            float tanh_y5 = tanh(fbk*(y(5) - b));
            
            // F(y) using current y values
            F(0) = g*(tanh_x-tanh_y1) + s1 - y(0);
            F(1) = g*(tanh_y1-tanh_y2) + s2 - y(1);
            F(2) = g*(tanh_y2-tanh_y3) + s3 - y(2);
            F(3) = g*(tanh_y3-tanh_y4) + s4 - y(3);
            F(4) = tanh_y5 - (s4 / (1.f - g_hp)) - y(4);
            
            // pre-compute re-used algebra (helper "functions")
            float help_x = 1.f - (tanh_x * tanh_x);
            float help_y1 = 1.f - (tanh_y1 * tanh_y1);
            float help_y2 = 1.f - (tanh_y2 * tanh_y2);
            float help_y3 = 1.f - (tanh_y3 * tanh_y3);
            float help_y4 = 1.f - (tanh_y4 * tanh_y4);
            float help_y5 = 1.f - (tanh_y5 * tanh_y5);
            float minus_g = -1.f * g;
            
            // Jacobian Matrix
            /*
            J = {
                {(minus_g*help_y1 - 1.0), 0.0, 0.0, (minus_g*r*help_x), (g*help_x)},
                {(g*help_y1), (minus_g*help_y2-1.0), 0.0, 0.0, 0.0},
                {0.0, (g*help_y2), (minus_g*help_y3 - 1.0), 0.0, 0.0},
                {0.0, 0.0, (g*help_y3), (minus_g*help_y4-1.0), 0.0},
                {0.0, 0.0, 0.0, (fbk*help_y5), -1.0}
            };
            */
            // only change these matrix elements (the rest are zeros)
            J(0,0) = minus_g*help_y1 - 1.f;
            J(0,3) = minus_g*r*help_x;
            J(0,4) = g*help_x;
            J(1,0) = g*help_y1;
            J(1,1) = minus_g*help_y2 - 1.f;
            J(2,1) = g*help_y2;
            J(2,2) = minus_g*help_y3 - 1.f;
            J(3,2) = g*help_y3;
            J(3,3) = minus_g*help_y4 - 1.f;
            J(4,3) = fbk*help_y5;
            J(4,4) = -1.f;
            
            // calculate next NR step
            y_est = y;
            delta_y = solve(J, (-1.f*F));
            y = delta_y + y;
            residue = y(3) = y_est(3);
            i++;
        }
        
        // update capacitor states
        s1 = 2.f * y(0) - s1;
        s2 = 2.f * y(1) - s2;
        s3 = 2.f * y(2) - s3;
        s4 = 2.f * y(3) - s4;
        s5 = 2.f * y(4) - s5;
        
        // save y
        //output = y(3);
        return y(3);
    }
};

MIN_EXTERNAL(ladder);

