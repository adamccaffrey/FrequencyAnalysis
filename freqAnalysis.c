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
#define     NORM 1.f / (3.f * (N/2.f))
#define		number_harmonics 20

static t_class *freqAnalysis_tilde_class;

typedef struct _freqAnalysis_tilde
{
	t_object		x_obj;

	float			x_signal_in, tone, real, prevPhase[N], maxBinValue, 
					maxBinCentreFrequency, dynamicNorm, pitchThresh;

	int 			harmonics, maxBinIndex, sampleCount, hopCount;

	t_float			sr, A[number_harmonics], inputLive[N], inputIR[N], analysisFreqs[N], analysisMags[N],
					liveFFTinput[N], irFFTinput[N], IFFToutput[N], output[N * OVERLAP];
	
	fftwf_complex	liveFFToutput[STEP], overtonesPol[STEP], overtonesCart[STEP], irFFToutput[STEP], IFFTinput[STEP];
	
	fftwf_plan      liveFFT, irFFT, IFFT;						

	t_inlet*		in2;
	t_outlet*		out;
	t_outlet*		out2;
	
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
float map(float x, float in_floor, float in_ceil, float out_floor, float out_ceil)
{
	return (x - in_floor) * (out_ceil - out_floor) / (in_ceil - in_floor) + out_floor;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static void freqAnalysis_tilde_harmonics(t_freqAnalysis *x, t_float f)
{
	f = f > 1 ? 1 : f < 1/number_harmonics ? 1/number_harmonics : f;

	x->harmonics = (int)((float)number_harmonics * f);
	x->harmonics = x->harmonics < 1 ? 1: x->harmonics > number_harmonics ? number_harmonics: x->harmonics;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static void freqAnalysis_tilde_tone(t_freqAnalysis *x, t_float f)
{
	f = f > 1 ? 1 : f < 0 ? 0 : f;
	x->tone = f;
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
static inline float profile(float fi, float bwi)
{
	float x = fi/bwi;
	x *= x;
	if (x>14.71280603) return 0.0;
    return (exp(-x)/bwi);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static void freqAnalysis_tilde_pitchThresh(t_freqAnalysis *x, t_float f)
{
	f = f > 1 ? 1 : f < 0 ? 0 : f;

	x->pitchThresh = f * N/3.0;
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
			if (x->maxBinValue > x->pitchThresh)
				x->maxBinIndex = i;
		}
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static void overtoneGen(t_freqAnalysis *x, float f, float bw, int nh)
{
	float bw_Hz = 0;
	float bwi = 0;
	float fi = 0;
	float hprofile = 0;

	for (nh = 1; nh < x->harmonics; nh++) {
		bw_Hz = (powf(2, bw / 1200.0) - 1.0) * f * (float)nh;
		bwi = bw_Hz / (2.0 * x->sr);
		fi = f * (float)nh / x->sr;
		for (int i = 0; i < STEP; i++) {
			hprofile = profile(((float)i / (float)N) - fi, bwi);
			x->overtonesPol[i][0] += hprofile * (float)x->A[nh];
		}
	}

	for (int i = 0; i < STEP; i++){
		x->overtonesPol[i][1] = (float)rand() / 2147483647.0 * 2.0 * M_PI;
		x->overtonesCart[i][0] = x->maxBinValue * NORM/10.0 * x->overtonesPol[i][0] * cosf(x->overtonesPol[i][1]);
		x->overtonesCart[i][1] = x->maxBinValue * NORM/10.0 * x->overtonesPol[i][0] * sinf(x->overtonesPol[i][1]);
	}
} 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static void freqAnalysis_tilde_bang(t_freqAnalysis *x)
{
	// for (int i = 0; i<STEP; i++)
	// {
	// 	// if (i % STEP/4 == 0)
	// 	// {
	// 	post("\nIndex: %d", i);
	// 	post("Overtone Re: %f", x->overtonesCart[i][0]);
	// 	post("Overtone Im: %f", x->overtonesCart[i][1]);
	// 	post("IFFTinput Re: %f", x->IFFTinput[i][0]);
	// 	post("IFFTinput Im: %f", x->IFFTinput[i][1]);
		
	// 	// }
	// }

	// post("MaxBinValue: %f", x->maxBinValue);
	// post("pitchThresh: %f", x->pitchThresh);
	post("Harmonics: %d", x->harmonics);
	// post("Max Frequency: %f", x->analysisFreqs[x->maxBinIndex]*x->sr/N);
	// post("binCentreFrequency: %f", x->maxBinCentreFrequency * x->sr/(2.0 * M_PI));
	// }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static t_int *freqAnalysis_tilde_perform(t_int *w)
{
	t_freqAnalysis* x = (t_freqAnalysis*)(w[1]);
	t_sample* in1 = (t_sample*)(w[2]);
	t_sample* in2 = (t_sample*)(w[3]);
	t_sample* out = (t_sample*)(w[4]);
	t_sample* out2 = (t_sample*)(w[5]);
	int blksize = (int)(w[6]);
	float bw = 50;// bandwidth of first harmonic in cents
	int nh = number_harmonics;// number of harmonics, nh < (sr/f)
	// float f = x->analysisFreqs[x->maxBinIndex]*x->sr/N;
	while(blksize--)
	{
		if(x->hopCount == 0)
		{
			//unwrap input
			for(int i = 0; i < N; i++) 
			{
				x->liveFFTinput[i] = (x->inputLive[(x->sampleCount+i)%N] * window((float)i/N));
			}
			//do FFT
			fftwf_execute(x->liveFFT);

			pitchDetect(x, x->liveFFToutput);

			overtoneGen(x, x->analysisFreqs[x->maxBinIndex]*x->sr/(float)N, bw, nh);

			for (int i=0; i<STEP; i++)
			{
				x->IFFTinput[i][0] = (1.0-x->tone/2) * x->liveFFToutput[i][0] + x->tone * x->overtonesCart[i][0];
				x->IFFTinput[i][1] = (1.0-x->tone/2) * x->liveFFToutput[i][1] + x->tone * x->overtonesCart[i][1];
			}
			//do Spectral Stuff  
			fftwf_execute(x->IFFT);
			//overlap-add result
			int i;
			int startix = 0;
			for(i = 0; i < N; i++) 
			{
				startix = (x->sampleCount + i) % N;
				if((N-i) <= HOP) x->output[startix] = 	(x->IFFToutput[i] * window((float)i/N));
				else 			 x->output[startix] += 	(x->IFFToutput[i] * window((float)i/N));
			}
		}

		//read/write input/output samples
		float liveInput = *in1;
		float irInput	= *in2;
		x->inputLive[x->sampleCount] = liveInput;
		x->inputIR[x->sampleCount] = irInput;

		// *out++ = (x->output[x->sampleCount] * NORM);
		*out++ = x->analysisFreqs[x->maxBinIndex]*x->sr/(float)N;

		//ignore this just incrementing unused i/o pointers
		++out2;
		++in1;
		++in2; 

		//increment and wrap counters
		x->sampleCount++;
		x->sampleCount %= N;
		x->hopCount = x->sampleCount % HOP;
	}
	return w+7;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static void *freqAnalysis_tilde_new()	
{
	t_freqAnalysis *x = (t_freqAnalysis *)pd_new(freqAnalysis_tilde_class);
	x->sr = sys_getsr();
	x->in2 = inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal);
	x->out = outlet_new(&x->x_obj, &s_signal);
	x->out2 = outlet_new(&x->x_obj, &s_signal);
	
	x->liveFFT = fftwf_plan_dft_r2c_1d(N, x->liveFFTinput, x->liveFFToutput, FFTW_MEASURE);		//real to complex fft plan (fft)
	x->irFFT = fftwf_plan_dft_r2c_1d(N, x->irFFTinput, x->irFFToutput, FFTW_MEASURE);		//real to complex fft plan (fft)
	x->IFFT = fftwf_plan_dft_c2r_1d(N, x->IFFTinput, x->IFFToutput, FFTW_MEASURE);		//complex to real fft plan (ifft)
	// x->IFFT = fftwf_plan_dft_c2r_1d(N, x->overtones, x->IFFToutput, FFTW_MEASURE);		//complex to real fft plan (ifft)
	
	x->harmonics = 2;

	x->dynamicNorm = NORM + (50.0 * NORM);

	x->A[0]=0.0;//A[0] is not used
    for (int i=1;i<number_harmonics;i++) x->A[i]=1.0/i;

	for (int i = 0; i < STEP; i++)
	{		
		x->liveFFToutput[i][0] = 0;
		x->liveFFToutput[i][1] = 0;
		x->irFFToutput[i][0] = 0;
		x->irFFToutput[i][1] = 0;
		x->IFFTinput[i][0] = 0;
		x->IFFTinput[i][1] = 0;
		x->overtonesPol[i][0] = 0;
		x->overtonesPol[i][1] = 0;
		x->overtonesCart[i][0] = 0;
		x->overtonesCart[i][1] = 0;
    }
	
	for(int i = 0; i < N; i++)
	{
		x->inputLive[i] = 0;
		x->inputIR[i] = 0;
		x->liveFFTinput[i] = 0;
		x->irFFTinput[i] = 0;
		x->IFFToutput[i] = 0;
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
	dsp_add(freqAnalysis_tilde_perform, 6, x, sp[0]->s_vec, sp[1]->s_vec, sp[2]->s_vec, sp[3]->s_vec, sp[0]->s_n);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static void freqAnalysis_tilde_free(t_freqAnalysis *x)
{
	inlet_free(x->in2);
	outlet_free(x->out);
	outlet_free(x->out2);	
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

	class_addmethod(freqAnalysis_tilde_class, (t_method)freqAnalysis_tilde_harmonics, gensym("harmonics"), A_FLOAT, 0);

	class_addmethod(freqAnalysis_tilde_class, (t_method)freqAnalysis_tilde_tone, gensym("tone"), A_FLOAT, 0);

	class_addmethod(freqAnalysis_tilde_class, (t_method)freqAnalysis_tilde_pitchThresh, gensym("pitchThresh"), A_FLOAT, 0);

}
//END//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
