#include "c74_min.h"
#include "../shared/signal_routing_objects.h"
#include <math.h>


class ladder : public object<ladder>, public sample_operator<4, 1>
{
private:
    // vector declaration
    float y[5] = {0.f, 0.f, 0.f, 0.f, 0.f};
    float y_est[5] = {0.f, 0.f, 0.f, 0.f, 0.f};
    float F[5] = {0.f, 0.f, 0.f, 0.f, 0.f};
    
    // jacobian elements
    float j00 = 0.f;
    float j03 = 0.f;
    float j04 = 0.f;
    float j10 = 0.f;
    float j11 = 0.f;
    float j21 = 0.f;
    float j22 = 0.f;
    float j32 = 0.f;
    float j33 = 0.f;
    float j43 = 0.f;
    float den = 0.f;
    
    // max NR error (stopping condition)
    float errorThresh = 0.000001;
        
    // accumulator states (capacitors)
    float s1 = 0.f;
    float s2 = 0.f;
    float s3 = 0.f;
    float s4 = 0.f;
    float s5 = 0.f;
    
    
public:
    // metadata
	MIN_DESCRIPTION {"Transistor ladder filter with feedback."};
	MIN_TAGS {"filter, ladder, lpf"};
	MIN_AUTHOR {"Peter Corbett (Osprey Instruments)"};
	MIN_RELATED {"lores~, svf~, onepole~"};

    // IO
	inlet<>  x {this, "(signal) Signal Input"};
	inlet<>  cutoff {this, "(signal) Cutoff Frequency"};
	inlet<>  r {this, "(signal) Resonance"};
    inlet<>  fbk {this, "(number) Feedback"};
	outlet<> output {this, "(signal) Output", "signal"};
    attribute<float> b { this, "bias", 0.65, range {-1.0, 1.0} };
    
    // more efficient tanh?
    float tanh2(float x)
    {
        float x2 = x*x;
        float sh = x*(1+x2*((1/6.f) + x2*(1/120.f)));
        return sh / sqrt(1+sh*sh);

    }

    // highpass filter
    float Fc_hp = 8.f;
    float g_hp = Fc_hp * M_PI / samplerate();
    float g_den = 2.f * g_hp + 1.f;
    
        
	// Call operator: process a single sample
	sample operator()(sample x, sample cutoff, sample r, sample fbk)
    {
        
        float g = cutoff * M_PI / samplerate();
        
        int i = 0;
        int residue = 100;
        while(abs(residue) > errorThresh && i < 50)
        {
            // pre-compute tanh functions
            float tanh_x = tanh(x - r*y[3] + y[4]);
            float tanh_y1 = tanh(y[0]);
            float tanh_y2 = tanh(y[1]);
            float tanh_y3 = tanh(y[2]);
            float tanh_y4 = tanh(y[3]);
            float tanh_y5 = tanh(fbk*(y[3] - b));
            
            // F(y) using current y values
            F[0] = g*(tanh_x-tanh_y1) + s1 - y[0];
            F[1] = g*(tanh_y1-tanh_y2) + s2 - y[1];
            F[2] = g*(tanh_y2-tanh_y3) + s3 - y[2];
            F[3] = g*(tanh_y3-tanh_y4) + s4 - y[3];
            F[4] = (tanh_y5 + s5) / g_den - y[4];
            
            // pre-compute re-used algebra (helper "functions")
            float help_x = 1.f - (tanh_x * tanh_x);
            float help_y1 = 1.f - (tanh_y1 * tanh_y1);
            float help_y2 = 1.f - (tanh_y2 * tanh_y2);
            float help_y3 = 1.f - (tanh_y3 * tanh_y3);
            float help_y4 = 1.f - (tanh_y4 * tanh_y4);
            float help_y5 = 1.f - (tanh_y5 * tanh_y5);
            float minus_g = -1.f * g;
            
            // Dynamic Jacobian Matrix Elements
            // only change these, the rest are zero (sparse)
            j00 = minus_g*help_y1 - 1.f;
            j03 = minus_g*r*help_x;
            j04 = g*help_x;
            j10 = g*help_y1;
            j11 = minus_g*help_y2 - 1.f;
            j21 = g*help_y2;
            j22 = minus_g*help_y3 - 1.f;
            j32 = g*help_y3;
            j33 = minus_g*help_y4 - 1.f;
            j43 = fbk*help_y5/g_den;

            
            // calculate next NR step
            std::copy(std::begin(y), std::end(y), std::begin(y_est));

            den = j00*j11*j22*j33 - j03*j10*j21*j32 - j04*j10*j21*j32*j43;
            y[0] = y[0] - (F[0]*j11*j22*j33 - F[1]*(j03*j21*j32 + j04*j21*j32*j43) + F[2]*(j03*j11*j32 + j04*j11*j32*j43) - F[3]*(j03*j11*j22 + j04*j11*j22*j43) + F[4]*j04*j11*j22*j33) / den;
            y[1] = y[1] + (F[0]*j10*j22*j33 - F[1]*j00*j22*j33 + F[2]*(j03*j10*j32 + j04*j10*j32*j43) - F[3]*(j03*j10*j22 + j04*j10*j22*j43) + F[4]*j04*j10*j22*j33) / den;
            y[2] = y[2] - (F[0]*j10*j21*j33 - F[1]*j00*j21*j33 + F[2]*j00*j11*j33 - F[3]*(j03*j10*j21 + j04*j10*j21*j43) + F[4]*j04*j10*j21*j33) / den;
            y[3] = y[3] + (F[0]*j10*j21*j32 - F[1]*j00*j21*j32 + F[2]*j00*j11*j32 - F[3]*j00*j11*j22 + F[4]*j04*j10*j21*j32) / den;
            y[4] = y[4] + (F[0]*j10*j21*j32*j43 - F[1]*j00*j21*j32*j43 + F[2]*j00*j11*j32*j43 - F[3]*j00*j11*j22*j43 + F[4]*(j00*j11*j22*j33 - j03*j10*j21*j32)) / den;
            
            residue = y[3] - y_est[3];
            i++;
        }
        
        // update capacitor states
        s1 = 2.f * y[0] - s1;
        s2 = 2.f * y[1] - s2;
        s3 = 2.f * y[2] - s3;
        s4 = 2.f * y[3] - s4;
        s5 = 2.f * (y[4] - tanh(fbk*(y[3] - b))) - s5;
        
        // save y
        return y[3];
        
	}
};

MIN_EXTERNAL(ladder);
