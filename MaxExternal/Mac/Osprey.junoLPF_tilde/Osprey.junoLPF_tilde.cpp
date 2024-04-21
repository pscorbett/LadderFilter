#include "c74_min.h"
#include "../shared/signal_routing_objects.h"
#include <math.h>


class junoLPF : public object<junoLPF>, public sample_operator<3, 4>
{
private:
    // vector declaration
    float y[4] = {0.f, 0.f, 0.f, 0.f};
    float y_est[4] = {0.f, 0.f, 0.f, 0.f};
    float F[4] = {0.f, 0.f, 0.f, 0.f};
    
    // jacobian elements
    float j00 = 0.f;
    float j03 = 0.f;
    float j10 = 0.f;
    float j11 = 0.f;
    float j21 = 0.f;
    float j22 = 0.f;
    float j32 = 0.f;
    float j33 = 0.f;
    float den = 0.f;
    
    // max NR error (stopping condition)
    float errorThresh = 0.000001;
        
    // accumulator states (capacitors)
    float s1 = 0.f;
    float s2 = 0.f;
    float s3 = 0.f;
    float s4 = 0.f;
    
    
public:
    // metadata
	MIN_DESCRIPTION {"OTA Cascade Filter"};
	MIN_TAGS {"cascade, Roland, Juno, lpf"};
	MIN_AUTHOR {"Peter Corbett (Osprey Instruments)"};
	MIN_RELATED {"lores~, svf~, onepole~"};

    // IO
	inlet<>  x {this, "(signal) Signal Input"};
	inlet<>  cutoff {this, "(signal) Cutoff Frequency"};
	inlet<>  r {this, "(signal) Resonance"};
	outlet<> out1p {this, "(signal) 1-Pole Output", "signal"};
    outlet<> out2p {this, "(signal) 2-Pole Output", "signal"};
    outlet<> out3p {this, "(signal) 3-Pole Output", "signal"};
    outlet<> out4p {this, "(signal) 4-Pole Output", "signal"};
    
    // more efficient tanh?
    float tanh2(float x)
    {
        float x2 = x*x;
        float sh = x*(1+x2*((1/6.f) + x2*(1/120.f)));
        return sh / sqrt(1+sh*sh);

    }

    
        
	// Call operator: process a single sample
	samples<4> operator()(sample x, sample cutoff, sample r)
    {
        
        float g = tan(cutoff * M_PI / samplerate());
        
        int i = 0;
        int residue = 100;
        while(abs(residue) > errorThresh && i < 50)
        {
            // pre-compute tanh functions
            float tanh_y1 = tanh(x - y[0] - r*y[3]);
            float tanh_y2 = tanh(y[0] - y[1]);
            float tanh_y3 = tanh(y[1] - y[2]);
            float tanh_y4 = tanh(y[2] - y[3]);
            
            // F(y) using current y values
            F[0] = g*tanh_y1 + s1 - y[0];
            F[1] = g*tanh_y2 + s2 - y[1];
            F[2] = g*tanh_y3 + s3 - y[2];
            F[3] = g*tanh_y4 + s4 - y[3];

            
            // pre-compute re-used algebra (helper "functions")
            float help_y1 = 1.f - (tanh_y1 * tanh_y1);
            float help_y2 = 1.f - (tanh_y2 * tanh_y2);
            float help_y3 = 1.f - (tanh_y3 * tanh_y3);
            float help_y4 = 1.f - (tanh_y4 * tanh_y4);
            float minus_g = -1.f * g;
            
            // Dynamic Jacobian Matrix Elements
            // only change these, the rest are zero (sparse)
            j00 = minus_g*help_y1 - 1.f;
            j03 = minus_g*r*help_y1;
            j10 = g*help_y2;
            j11 = minus_g*help_y2 - 1.f;
            j21 = g*help_y3;
            j22 = minus_g*help_y3 - 1.f;
            j32 = g*help_y4;
            j33 = minus_g*help_y4 - 1.f;
            den = j00*j11*j22*j33 - j03*j10*j21*j32;

            
            // calculate next NR step
            std::copy(std::begin(y), std::end(y), std::begin(y_est));

            y[0] = y[0] + (F[1]*j03*j21*j32 - F[0]*j11*j22*j33 - F[2]*j03*j11*j32 + F[3]*j03*j11*j22) / den;
            y[1] = y[1] + (F[0]*j10*j22*j33 - F[1]*j00*j22*j33 + F[2]*j03*j10*j32 - F[3]*j03*j10*j22) / den;
            y[2] = y[2] + (F[1]*j00*j21*j33 - F[0]*j10*j21*j33 - F[2]*j00*j11*j33 + F[3]*j03*j10*j21) / den;
            y[3] = y[3] + (F[0]*j10*j21*j32 - F[1]*j00*j21*j32 + F[2]*j00*j11*j32 - F[3]*j00*j11*j22) / den;

            
            residue = y[3] - y_est[3];
            i++;
        }
        
        // update capacitor states
        s1 = 2.f * y[0] - s1;
        s2 = 2.f * y[1] - s2;
        s3 = 2.f * y[2] - s3;
        s4 = 2.f * y[3] - s4;

        
        // save y
        return {{y[0], y[1], y[2], y[3]}};
        
	}
};

MIN_EXTERNAL(junoLPF);
