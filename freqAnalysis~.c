//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// freqAnalysis~																																																	//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include	"m_pd.h"
#include	<fftw3.h>
#include	<stdbool.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>
#include	"windows.h"

#define     N 2048		//fft
#define 	STEP ((N/2)+1)
#define 	OVERLAP 4
#define		HOP (N/OVERLAP)
// #define     NORM 2.f/(3.f*N)
#define     NORM 1.f/(100.f*N)

static t_class *freqAnalysis_tilde_class;

typedef struct _freqAnalysis_tilde
{
	t_object		x_obj;

	float			x_signal_in, real, prevPhase[N], maxBinValue, env, attack, decay, hold, filterFreq, filterQ, 
					maxBinCentreFrequency, dynamicNorm, thresh, outputFrequency;

	int 			maxBinIndex, sampleCount, hopCount;

	t_float			sr, inputLive[N], analysisFreqs[N], analysisMags[N],
					liveFFTinput[N], output[N * OVERLAP];
	
	fftwf_complex	liveFFToutput[STEP];
	
	fftwf_plan      liveFFT;						

	t_inlet*		in2;
	t_outlet*		out;
	t_outlet*		out2;
	t_outlet*		out3;

	} t_freqAnalysis;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static inline float window(t_float position)
{
    float fpos = position * N;
    uint16_t ipos = fpos;
    float frac = fpos - ipos;
    return (hann[ipos]+(frac*(hann[ipos]-hann[ipos+1])));
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
float peak(fftwf_complex *input)
{
    float mag, peak;
    peak = 0;
    for(int i = 0; i < STEP; i++)  
    {
        mag = sqrtf(input[i][0] * input[i][0]) + (input[i][1] * input[i][1]);
        if(mag > peak) peak = mag;
    }
    return peak;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
float rms(fftwf_complex *input)
{
    float totalMag, rms;
    totalMag = 0;
    for(int i = 0; i < STEP; i++) totalMag += sqrtf(input[i][0] * input[i][0]) + (input[i][1] * input[i][1]);
    rms = totalMag / STEP;
    return rms;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static inline float map(float x, float in_floor, float in_ceil, float out_floor, float out_ceil)
{
	return (x - in_floor) * (out_ceil - out_floor) / (in_ceil - in_floor) + out_floor;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static inline float modAlt(float x, float y)
{
    return x - (y * (float)((int)(x/y)));
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static inline float wrapPhase(float phaseIn)
{
    if (phaseIn >= 0)
        return modAlt(phaseIn + M_PI, 2.0 * M_PI) - M_PI;
    else
        return modAlt(phaseIn - M_PI, -2.0 * M_PI) + M_PI;	
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static void freqAnalysis_tilde_thresh(t_freqAnalysis *x, t_float f)
{
	f = f > 1 ? 1 : f < 0 ? 0 : f;
	x->thresh = f * N/3.0;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static void freqAnalysis_tilde_hold(t_freqAnalysis *x, t_float f)
{
	f = f > 1 ? 1 : f < 0 ? 0 : f;
	x->hold = f;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static void freqAnalysis_tilde_attack(t_freqAnalysis *x, t_float f)
{
	f = f > 1 ? 1 : f < 0 ? 0 : 1.f - f;
	x->attack = 1.f - (f * f * f * 0.999999f);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static void freqAnalysis_tilde_decay(t_freqAnalysis *x, t_float f)
{
	f = f > 1 ? 1 : f < 0 ? 0 : (1.f - f);
	x->decay = 1.f - (f * f * f * 0.1);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static void freqAnalysis_tilde_filterFreq(t_freqAnalysis *x, t_float f)
{
	f = f > 1 ? 1 : f < 0 ? 0 : f;
	x->filterFreq = f;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static void freqAnalysis_tilde_filterQ(t_freqAnalysis *x, t_float f)
{
	f = f > 1 ? 1 : f < 0 ? 0 : f;
	x->filterQ = f;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static void pitchDetect(t_freqAnalysis *x, fftwf_complex input[STEP])
{
	float real = 0;
	float imag = 0;
	float mag = 0;
	float phase = 0;
	float phaseDiff = 0;
	float binDeviation = 0;
	float binCentreFrequency = 0;
	float maxBinValue = 0;
	
	for (int i = 0; i < STEP; i++)
	{
		real = input[i][0];
		imag = input[i][1];
		mag = sqrtf((real * real) + (imag * imag));

		if (real < 0 && imag < 0)
		{ phase = atan2f(imag - M_PI, real); }
		if (real < 0 && imag > 0)
		{ phase = atan2f(imag + M_PI, real); }
		else{ phase = atan2f(imag, real); }

		phaseDiff = phase - x->prevPhase[i]; 
		binCentreFrequency = 2.0 * M_PI * (float)i / (float)N;
		phaseDiff = wrapPhase(phaseDiff - binCentreFrequency * HOP);
		binDeviation = phaseDiff * (float)N / (float)HOP  / (2.0 * M_PI);
		x->analysisFreqs[i] = (float)i + binDeviation;
		x->analysisMags[i] = mag;
		x->prevPhase[i] = phase;
		
		if (mag > maxBinValue){
			maxBinValue = mag;
			x->maxBinValue = maxBinValue;
			x->maxBinCentreFrequency = binCentreFrequency;
			if (x->maxBinValue > x->thresh)
				x->maxBinIndex = i;
		}
	}
	x->outputFrequency = x->analysisFreqs[x->maxBinIndex]*x->sr/(float)N;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static void freqAnalysis_tilde_bang(t_freqAnalysis *x)
{
	post("Attack: %f", x->attack);
	post("Decay: %f", x->decay);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static t_int *freqAnalysis_tilde_perform(t_int *w)
{
	t_freqAnalysis* x = (t_freqAnalysis*)(w[1]);
	t_sample* in1 = (t_sample*)(w[2]);
	t_sample* in2 = (t_sample*)(w[3]);
	t_sample* out = (t_sample*)(w[4]);
	t_sample* out2 = (t_sample*)(w[5]);
	t_sample* out3 = (t_sample*)(w[6]);

	int blksize = (int)(w[7]);
	float input_peak = 0.f;
	float input_rms = 0.f;
	float delta_peak = 0.f;
	float delta_rms = 0.f;

	while(blksize--)
	{
		if(x->hopCount == 0)
		{
			for(int i = 0; i < N; i++) 
			{
				x->liveFFTinput[i] = (x->inputLive[(x->sampleCount+i)%N] * window((float)i/N));
			}
			fftwf_execute(x->liveFFT);
			pitchDetect(x, x->liveFFToutput);

			input_peak = peak(x->liveFFToutput);
			input_rms = rms(x->liveFFToutput);

			delta_peak = input_peak * NORM - x->env;
			delta_rms = input_rms - x->env;

			if (input_peak * NORM >= x->env){
				x->env = x->attack * (delta_rms + x->env);
			}
			else {
				x->env = x->decay * (x->env - delta_peak) + delta_peak;
			}
			x->env = x->env < 0 ? 0 : x->env > 1 ? 1: x->env;

		}

		float liveInput = *in1;
		x->inputLive[x->sampleCount] = liveInput;

		*out++ = x->outputFrequency;
		*out2++ = x->env;
		*out3++ = input_rms * NORM;
		++in1;
		++in2; 

		//increment and wrap counters
		x->sampleCount++;
		x->sampleCount %= N;
		x->hopCount = x->sampleCount % HOP;
	}
	return w+8;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static void *freqAnalysis_tilde_new()	
{
	t_freqAnalysis *x = (t_freqAnalysis *)pd_new(freqAnalysis_tilde_class);
	x->sr = sys_getsr();
	x->in2 = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
	x->out = outlet_new(&x->x_obj, &s_signal);
	x->out2 = outlet_new(&x->x_obj, &s_signal);
	x->out3 = outlet_new(&x->x_obj, &s_signal);

	x->liveFFT = fftwf_plan_dft_r2c_1d(N, x->liveFFTinput, x->liveFFToutput, FFTW_MEASURE);		//real to complex fft plan (fft)

	x->outputFrequency = 0.f;
	x->attack = 1.f;
	x->decay = 1.f;
	x->hold = 0.f;
	x->filterFreq = 0.f;
	x->filterQ = 0.f;
	x->thresh = 0.f;
	x->env = 0;


	x->dynamicNorm = NORM + (50.0 * NORM);

	for (int i = 0; i < STEP; i++)
	{		
		x->liveFFToutput[i][0] = 0;
		x->liveFFToutput[i][1] = 0;
    }
	
	for(int i = 0; i < N; i++)
	{
		x->inputLive[i] = 0;
		x->liveFFTinput[i] = 0;
		x->output[i] = 0;
		x->prevPhase[i] = 0;
		x->analysisFreqs[i] = 0;
		x->analysisMags[i] = 0;
	}
 
	return (void *)x;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static void freqAnalysis_tilde_dsp(t_freqAnalysis *x, t_signal **sp)
{
	dsp_add(freqAnalysis_tilde_perform, 7, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, sp[4]->s_vec, sp[0]->s_n);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static void freqAnalysis_tilde_free(t_freqAnalysis *x)
{
	inlet_free(x->in2);
	outlet_free(x->out);
	outlet_free(x->out2);	
	outlet_free(x->out3);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void freqAnalysis_tilde_setup(void)
{
	freqAnalysis_tilde_class = class_new(gensym("freqAnalysis~"),
										(t_newmethod)freqAnalysis_tilde_new,
										(t_method)freqAnalysis_tilde_free,
										sizeof(t_freqAnalysis),
										CLASS_DEFAULT,
										0);
	CLASS_MAINSIGNALIN(freqAnalysis_tilde_class, t_freqAnalysis, x_signal_in);

	class_addbang(freqAnalysis_tilde_class, (t_method)freqAnalysis_tilde_bang);

	class_addmethod(freqAnalysis_tilde_class, (t_method)freqAnalysis_tilde_dsp, gensym("dsp"), A_CANT, 0);

	class_addmethod(freqAnalysis_tilde_class, (t_method)freqAnalysis_tilde_thresh, gensym("thresh"), A_FLOAT, 0);

	class_addmethod(freqAnalysis_tilde_class, (t_method)freqAnalysis_tilde_hold, gensym("hold"), A_FLOAT, 0);

	class_addmethod(freqAnalysis_tilde_class, (t_method)freqAnalysis_tilde_attack, gensym("attack"), A_FLOAT, 0);

	class_addmethod(freqAnalysis_tilde_class, (t_method)freqAnalysis_tilde_decay, gensym("decay"), A_FLOAT, 0);

	class_addmethod(freqAnalysis_tilde_class, (t_method)freqAnalysis_tilde_filterFreq, gensym("filterFreq"), A_FLOAT, 0);
	
	class_addmethod(freqAnalysis_tilde_class, (t_method)freqAnalysis_tilde_filterQ, gensym("filterQ"), A_FLOAT, 0);

}
//END//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
