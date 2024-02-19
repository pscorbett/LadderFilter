#include "c74_min.h"
#include "../shared/signal_routing_objects.h"
#include <math.h>


class VCA_Filter : public object<VCA_Filter>, public sample_operator<3, 1>
{
private:
    // vector declaration
    float y[2] = {0.f, 0.f};
    float y_est[2] = {0.f, 0.f};
    float F[2] = {0.f, 0.f};
    
    // jacobian elements
    float j00 = 0.f;
    float j01 = 0.f;
    float j10 = 0.f;
    float j11 = 0.f;
    float den = 0.f;
    
    // max NR error (stopping condition)
    float errorThresh = 0.000001;
        
    // accumulator states (capacitors)
    float s1 = 0.f;
    float s2 = 0.f;

    
    
public:
    // metadata
	MIN_DESCRIPTION {"2-Pole Resonant VCA Filter"};
	MIN_TAGS {"filter, ladder, lpf"};
	MIN_AUTHOR {"Peter Corbett (Osprey Instruments)"};
	MIN_RELATED {"lores~, svf~, onepole~"};

    // IO
	inlet<>  x {this, "(signal) Signal Input"};
	inlet<>  cutoff {this, "(signal) Cutoff Frequency"};
	inlet<>  r {this, "(signal) Resonance"};
	outlet<> output {this, "(signal) Output", "signal"};
    attribute<float> b { this, "bias", 0.0, range {-1.0, 1.0} };
    
    // more efficient tanh?
    float tanh2(float x)
    {
        float x2 = x*x;
        float sh = x*(1+x2*((1/6.f) + x2*(1/120.f)));
        return sh / sqrt(1+sh*sh);

    }

    
        
	// Call operator: process a single sample
	sample operator()(sample x, sample cutoff, sample r)
    {
        
        float g = tan(cutoff * M_PI / samplerate());
        
        int i = 0;
        int residue = 100;
        while(abs(residue) > errorThresh && i < 50)
        {

            
            // F(y) using current y values
            F[0] = g*(x - y[0] - r*y[1]) + s1 - y[0];
            F[1] = g*(y[0] - y[1]) + s2 - y[1];
            
            // pre-compute re-used algebra (helper "functions")
            float jhelp = -1.f * (g + 1.f);
            
            // Dynamic Jacobian Matrix Elements
            // only change these, the rest are zero (sparse)
            j00 = jhelp;
            j01 = -1.f*g*r;
            j10 = g;
            j11 = jhelp;
            den = j00*j11 - j01*j10;


            
            // calculate next NR step
            std::copy(std::begin(y), std::end(y), std::begin(y_est));

            y[0] = y[0] + (F[1]*j01 - F[0]*j11) / den;
            y[1] = y[1] + (F[0]*j10 - F[1]*j00) / den;
            
            residue = y[1] - y_est[1];
            i++;
        }
        
        // update capacitor states
        s1 = 2.f * y[0] - s1;
        s2 = 2.f * y[1] - s2;

        
        // save y with OTA bias and nonlinearity
        return tanh(y[1] + b);
        
	}
};

MIN_EXTERNAL(VCA_Filter);
