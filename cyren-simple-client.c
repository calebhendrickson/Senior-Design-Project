/** @file simple_client.c
 *
 * @brief This simple client demonstrates the most basic features of JACK
 * as they would be used by many applications.
 */

#include <errno.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <pthread.h>
#include <strings.h>
#include <string.h>
#include <termios.h>
#include <math.h>
#include <time.h>
#include <sys/types.h>

#include <jack/jack.h>

#define M_PI 3.14159265358979323846

jack_port_t *input_port;
jack_port_t *output_port;
jack_client_t *client;

// tremolo external variables
float t = 0.0;

//user params
static float tremolo_rate = 0.5;  // (range: 0.0 to 1.0)
static float tremolo_depth = 0.5; // (range: 0.0 to 1.0)

// hard clipper
double max_hard_clipper = 0;

// overdrive
double max_overdrive = 0;

//chorus vars
int initChorus = 0;
float* chorusData; // Our own circular buffer of samples
int chorusBufLength = round(192000 * 0.02); // Length of our delay buffer in samples (Samplerate * amount of delay)

//delay vars
int initDelay = 0;
double* delayData;  // Our own circular buffer of samples   
int delayBufLength = round(192000 * 0.1); // Length of our delay buffer in samples (Samplerate * amount of delay) 

//flanger vars
int initFlanger = 0;
float* flangerData; // Our own circular buffer of samples
int flangerBufLength = round(0.003 * 192000); // Length of our delay buffer in samples (Samplerate * amount of delay) 


char *currentFilter;

typedef struct keyQueue {
    struct keyQueue *next;
    char key;
} keyQueue_t;

typedef struct ThreadInfo {
    pthread_t tid;           /* thread id */
    pthread_mutex_t kqmutex; /* protects key queue from race condition between threads */
    keyQueue_t kqhead;       /* input keys queued to this thread */
    char *keys;              /* keys this thread responds to */
    char *name;              /* name of this thread */
} threadInfo_t;

static struct termios origtc, newtc;

threadInfo_t threads[] = { 
	{ 0, PTHREAD_MUTEX_INITIALIZER, { NULL, '\0' }, "0", "No Filter" },
    { 0, PTHREAD_MUTEX_INITIALIZER, { NULL, '\0' }, "1", "Chorus" },
    { 0, PTHREAD_MUTEX_INITIALIZER, { NULL, '\0' }, "2", "Tremolo" },
	{ 0, PTHREAD_MUTEX_INITIALIZER, { NULL, '\0' }, "3", "Low Pass" },
    { 0, PTHREAD_MUTEX_INITIALIZER, { NULL, '\0' }, "4", "High Pass" },
    { 0, PTHREAD_MUTEX_INITIALIZER, { NULL, '\0' }, "5", "Flanger" },
    { 0, PTHREAD_MUTEX_INITIALIZER, { NULL, '\0' }, "6", "Hard Clipper" },
    { 0, PTHREAD_MUTEX_INITIALIZER, { NULL, '\0' }, "7", "Soft Clipper" },
    { 0, PTHREAD_MUTEX_INITIALIZER, { NULL, '\0' }, "8", "Fuzz" },
    { 0, PTHREAD_MUTEX_INITIALIZER, { NULL, '\0' }, "9", "Overdrive" },
    { 0, PTHREAD_MUTEX_INITIALIZER, { NULL, '\0' }, "-", "Delay" }
};

// TODO: Need to take out mux maybe? Don't worry about fighting for resources?

void *service(void *arg) {
    char key;
    threadInfo_t *t = &threads[(int)arg];    // get pointer to thread
    for(;;) {
        pthread_mutex_lock(&t->kqmutex);     // lock other threads out while we tamper 
        key = '\0';                          // initialize key to NULL
        if (t->kqhead.next != NULL) {        // Anything queued up for us?
            keyQueue_t *kq = t->kqhead.next; // if so get ptr to key pkt
            key = kq->key;                   // fetch key from pkt
            t->kqhead.next = kq->next;       // Point to next key in queue (or NULL if no more queued up).
            free(kq);
			if(strcmp(currentFilter, t->name) != 0){
		    	currentFilter = t->name;
				printf("'%c' was pressed.\nCurrent Filter: %s \n", key, currentFilter);
			}
        }  
        pthread_mutex_unlock(&t->kqmutex);   // unlock key queue

        // if (key != '\0') {                   // if we got a key, log it
        //     printf("'%c' was pressed.\nCurrent Filter: %s \n", key, currentFilter);
        // }

        // ⇓ usleep() probably more practical as 1-sec too long for most cases
        sleep(1);                            // sleep so we don't loop too fast eating CPU
    }
    return NULL;
}


/**
 * The process callback for this JACK application is called in a
 * special realtime thread once for each audio cycle.
 *
 * This client does nothing more than copy data from its input
 * port to its output port. It will exit when stopped by 
 * the user (e.g. using Ctrl-C on a unix-ish operating system)
 */
int
process (jack_nframes_t nframes, void *arg)
{
	jack_default_audio_sample_t *in, *out;
	
	in = jack_port_get_buffer (input_port, nframes);

	if(strcmp(currentFilter, "No Filter") == 0) {
	    //printf("%d", (int)nframes);
		out = jack_port_get_buffer (output_port, nframes);
		memcpy (out, in, sizeof (jack_default_audio_sample_t) * nframes);
	}
	else if(strcmp(currentFilter, "Chorus") == 0) {
        	   
    	jack_default_audio_sample_t *chorus;
		chorus = (jack_default_audio_sample_t *)malloc(sizeof(jack_default_audio_sample_t ) * nframes);
		memcpy (chorus, in, sizeof (jack_default_audio_sample_t) * nframes);
    	
    	// Variables whose values are set externally:
    	int numSamples = (int)nframes;     // How many audio samples to process
    	int dpw = 0;            // Write pointer into the delay buffer
    	float ph = 0.0;           // Current LFO phase, always between 0-1
    	float inverseSampleRate = 1 / 192000; // 1/f_s, where f_s = sample rate
    	
    	if(initChorus == 0){
    	    for(int i = 0; i < chorusBufLength; i++){
    	        *(chorusData + i) = 0;
    	    }
    	    initChorus = 1;
    	}
    									  
    	// User-adjustable effect parameters:
    	float frequency_ = 2;   // Frequency of the LFO //                      (range: 0 to 20)
    	float sweepWidth_ = 0.002 * 192000;  // Width of the LFO in samples     (range: unsure)
    	float depth_ = 0.5;       // Amount of delayed signal mixed with original (range: 0-1)
    	float feedback_ = 0.5;    // Amount of feedback                         (range: >= 0, < 1)    
    	
        for (int i = 0; i < numSamples; ++i)
        {
        	const float in = *(chorus + i);
        	float interpolatedSample = 0.0;
        	// Recalculate the read pointer position with respect to 
        	// the write pointer.
        	float currentDelay = sweepWidth_ * (0.5f + 0.5f * sinf(2.0 * M_PI * ph));
        	// Subtract 3 samples to the delay pointer to make sure 
        	// we have enough previous samples to interpolate with
        	float dpr = fmodf((float)dpw - (float)(currentDelay * 192000) + (float)chorusBufLength - 3.0, (float)chorusBufLength);
        	// Use linear interpolation to read a fractional index 
        	// into the buffer.
        	if(dpr < 0){
        	    float fraction = dpr - floorf(dpr);
        		float previousSample = floorf(dpr);
        		int nextSample = ((int)previousSample + 1) % chorusBufLength;
        		interpolatedSample = (fraction * *(chorusData + ((chorusBufLength + nextSample) % chorusBufLength))) + ((1.0f - fraction) * (*(chorusData + (chorusBufLength + (int)previousSample) % chorusBufLength)));
        	    *(chorusData + dpw) = in +  (interpolatedSample * feedback_);
        		// Increment the write pointer at a constant rate.
        		if (++dpw >= chorusBufLength){
        			dpw = 0;
        		}
        		// Store the output in the buffer, replacing the input
        		*(chorus + i) = in + depth_ * interpolatedSample;
        		// Update the LFO phase, keeping it in the range 0-1
        		ph += frequency_ * inverseSampleRate;
        		if (ph >= 1.0){
        			ph -= 1.0;
        		}
        	}
        	else{
        		// Use linear interpolation to read a fractional index 
        		// into the buffer.
        		float fraction = dpr - floorf(dpr);
        		int previousSample = (int)floorf(dpr);
        		int nextSample = (previousSample + 1) % chorusBufLength;
        		interpolatedSample = fraction * *(chorusData + nextSample) + (1.0f - fraction) * (*(chorusData + previousSample));
        		// Store the current information in the delay buffer. 
        		// With feedback, what we read is included in what gets 
        		// stored in the buffer, otherwise it’s just a simple 
        		// delay line of the input signal.
        		*(chorusData + dpw) = in + (interpolatedSample * feedback_);
        		// Increment the write pointer at a constant rate.
        		if (++dpw >= chorusBufLength){
        			dpw = 0;
        		}
        		// Store the output in the buffer, replacing the input
        		*(chorus + i) = in + depth_ * interpolatedSample;
        		// Update the LFO phase, keeping it in the range 0-1
        		ph += frequency_ * inverseSampleRate;
        		if (ph >= 1.0){
        			ph -= 1.0;
        		}
        	}
    	}
    	
    	out = jack_port_get_buffer (output_port, nframes);
		memcpy (out, chorus, sizeof (jack_default_audio_sample_t) * nframes);
		free(chorus);
		
    }
	else if(strcmp(currentFilter, "Tremolo") == 0) {

		jack_default_audio_sample_t *tremolo;
		tremolo = (jack_default_audio_sample_t *)malloc(sizeof(jack_default_audio_sample_t ) * nframes);
		memcpy (tremolo, in, sizeof (jack_default_audio_sample_t) * nframes);

        // no adjustable params

		for (int i = 0; i < (int)nframes; i++) {
			float trem_factor = (float)(1.0 - (tremolo_depth * (0.5 * sinf(t) + 0.5)));
			t += (float)(tremolo_rate * 0.002);

			if (t > 6.28318531) {
				t -= (float)6.28318531;
			}
			*(tremolo + i) = *(tremolo + i) * trem_factor;
		}

		out = jack_port_get_buffer (output_port, nframes);
		memcpy (out, tremolo, sizeof (jack_default_audio_sample_t) * nframes);
		free(tremolo);
	}
	else if(strcmp(currentFilter, "Low Pass") == 0) {
	    
	    jack_default_audio_sample_t *lowPass;
		lowPass = (jack_default_audio_sample_t *)malloc(sizeof(jack_default_audio_sample_t ) * nframes);
		memcpy (lowPass, in, sizeof (jack_default_audio_sample_t) * nframes);
		
		jack_default_audio_sample_t *lp_output;
		lp_output = (jack_default_audio_sample_t *)malloc(sizeof(jack_default_audio_sample_t ) * nframes);
		
        // param
        double d = 0.91; // seting the cutoff of the low pass (range: 0.9 to unsure)

    	for (int i = 0; i < (int)nframes; i++) {
    		if (i == 0) {
    			*(lp_output + i) = *(lowPass + i);
    		}
    		else {
    			*(lp_output + i) = ((1.0 - d) * (*(lowPass + i))) + (d * (*(lp_output + i - 1)));
    			// *(lp_output + i) = (*(lowPass + i) + (*(lowPass + i - 1)))/2;
    		}
    	}

		out = jack_port_get_buffer (output_port, nframes);
		memcpy (out, lp_output, sizeof (jack_default_audio_sample_t) * nframes);
		free(lp_output);
		free(lowPass);
		
	}
	else if(strcmp(currentFilter, "High Pass") == 0) {
	   
        jack_default_audio_sample_t *highPass;
		highPass = (jack_default_audio_sample_t *)malloc(sizeof(jack_default_audio_sample_t ) * nframes);
		memcpy (highPass, in, sizeof (jack_default_audio_sample_t) * nframes);
		
		jack_default_audio_sample_t *hp_output;
		hp_output = (jack_default_audio_sample_t *)malloc(sizeof(jack_default_audio_sample_t ) * nframes);
		
		// no adjustable params
		
    	for (int i = 0; i < (int)nframes; i++) {
    		if (i == 0) {
    			*(hp_output + i) = *(highPass + i);
    		}
    		else {
    			// *(hp_output + i) = (*(highPass + i) - (*(highPass + i - 1)))/2
    			*(hp_output + i) = (*(highPass + i) - (*(highPass + i - 1))) * 16;
    		}
    	}
    
		out = jack_port_get_buffer (output_port, nframes);
		memcpy (out, hp_output, sizeof (jack_default_audio_sample_t) * nframes);
		free(hp_output);
		free(highPass);
		
	}
	else if(strcmp(currentFilter, "Flanger") == 0) {
	    
	    jack_default_audio_sample_t *flanger;
		flanger = (jack_default_audio_sample_t *)malloc(sizeof(jack_default_audio_sample_t ) * nframes);
		memcpy (flanger, in, sizeof (jack_default_audio_sample_t) * nframes);
    	
    	// Variables whose values are set externally:
    	int numSamples = (int)nframes;     // How many audio samples to process
    	int dpw = 0;            // Write pointer into the delay buffer
    	float ph = 0.0;           // Current LFO phase, always between 0-1
    	float inverseSampleRate = 1 / 192000; // 1/f_s, where f_s = sample rate
    	
    	if(initFlanger == 0){
    	    for(int i = 0; i < flangerBufLength; i++){
    	        *(flangerData + i) = 0;
    	        // if this doesn't work try *(flangerData + i) = *(flanger + i);
    	    }
    	    initFlanger = 1;
    	}
    									  
    	// User-adjustable effect parameters:
    	float frequency_ = 2;   // Frequency of the LFO                         (range: 0 to 20hz)
    	float sweepWidth_ = 0.002 * 192000;  // Width of the LFO in samples     (range: unsure)
    	float depth_ = 0.5;       // Amount of delayed signal mixed with original (range: 0 to 1)
    	float feedback_ = 0.5;    // Amount of feedback                         (range: >= 0, < 1)    
    	
    	for (int i = 0; i < numSamples; ++i)
    	{
    		const float in = *(flanger + i);
    		float interpolatedSample = 0.0;
    		// Recalculate the read pointer position with respect to 
    		// the write pointer.
    		float currentDelay = sweepWidth_ * (0.5f + 0.5f * sinf(2.0 * M_PI * ph));
    		// Subtract 3 samples to the delay pointer to make sure 
    		// we have enough previous samples to interpolate with
    		float dpr = fmodf((float)dpw - (float)(currentDelay * 192000) + (float)flangerBufLength - 3.0, (float)flangerBufLength);
    		// Use linear interpolation to read a fractional index 
    		// into the buffer.
    		if(dpr < 0){
    		        float fraction = dpr - floorf(dpr);
            		float previousSample = floorf(dpr);
            		//int previousSample = (int)floorf(dpr);
            		int nextSample = ((int)previousSample + 1) % flangerBufLength;
            		interpolatedSample = (fraction * *(flangerData + ((flangerBufLength + nextSample) % chorusBufLength))) + ((1.0f - fraction) * (*(flangerData + ((flangerBufLength + (int)previousSample) % flangerBufLength))));
        		    *(flangerData + dpw) = in + (interpolatedSample * feedback_);
            		// Increment the write pointer at a constant rate.
            		if (++dpw >= flangerBufLength){
            			dpw = 0;
            		}
            		// Store the output in the buffer, replacing the input
            		*(flanger + i) = in + depth_ * interpolatedSample;
            		// Update the LFO phase, keeping it in the range 0-1
            		ph += frequency_ * inverseSampleRate;
            		if (ph >= 1.0){
            			ph -= 1.0;
            		}
        		}
        		else{
            		// Use linear interpolation to read a fractional index 
            		// into the buffer.
            		float fraction = dpr - floorf(dpr);
            		int previousSample = (int)floorf(dpr);
            		//int previousSample = (int)floorf(dpr);
            		int nextSample = (previousSample + 1) % flangerBufLength;
            		interpolatedSample = fraction * *(flangerData + nextSample) + (1.0f - fraction) * (*(flangerData + previousSample));
            		// Store the current information in the delay buffer. 
            		// With feedback, what we read is included in what gets 
            		// stored in the buffer, otherwise it’s just a simple 
            		// delay line of the input signal.
            		*(flangerData + dpw) = in + (interpolatedSample * feedback_);
            		// Increment the write pointer at a constant rate.
            		if (++dpw >= flangerBufLength){
            			dpw = 0;
            		}
            		// Store the output in the buffer, replacing the input
            		*(flanger + i) = in + depth_ * interpolatedSample;
            		// Update the LFO phase, keeping it in the range 0-1
            		ph += frequency_ * inverseSampleRate;
            		if (ph >= 1.0){
            			ph -= 1.0;
            		}
        		}
    	}
    	
    	out = jack_port_get_buffer (output_port, nframes);
		memcpy (out, flanger, sizeof (jack_default_audio_sample_t) * nframes);
		free(flanger);
	}
	else if(strcmp(currentFilter, "Hard Clipper") == 0) {
	    
	    jack_default_audio_sample_t *hard_clipper;
		hard_clipper = (jack_default_audio_sample_t *)malloc(sizeof(jack_default_audio_sample_t ) * nframes);
		memcpy (hard_clipper, in, sizeof (jack_default_audio_sample_t) * nframes);
	    
	    // finding max amplitude
	    for (int i = 0; i < (int)nframes; i++) {
    		if ((*(hard_clipper + i)) > max_hard_clipper) {
    			max_hard_clipper = (*(hard_clipper + i));
    		}
    	}
    	
    	// params
    	double hard_clipper_gain = 1.2;  // gain of the hard clipper (range: 1.0 to 1.5)
    	
    	
    	for (int i = 0; i < (int)nframes; i++) {
    		if (((*(hard_clipper + i)) * hard_clipper_gain) <= (-1 * max_hard_clipper)) {
    			*(hard_clipper + i) = -1 * max_hard_clipper;
    		}
    		else if (((*(hard_clipper + i)) * hard_clipper_gain) < max_hard_clipper) {
    			*(hard_clipper + i) = *(hard_clipper + i);
    		}
    		else {
    			*(hard_clipper + i) = max_hard_clipper;
    		}
    	}
    	
    	out = jack_port_get_buffer (output_port, nframes);
		memcpy (out, hard_clipper, sizeof (jack_default_audio_sample_t) * nframes);
		free(hard_clipper);
    	
	}
	else if(strcmp(currentFilter, "Soft Clipper") == 0) {
	    
	    // variables set externally
	    int sgn = 1;
        double e = 2.7182818;
        double temp = 0;
        
        //params
        double soft_clipper_gain = 3; // gain of the soft clipper (range: unsure)
	    
	    jack_default_audio_sample_t *soft_clipper;
		soft_clipper = (jack_default_audio_sample_t *)malloc(sizeof(jack_default_audio_sample_t ) * nframes);
		memcpy (soft_clipper, in, sizeof (jack_default_audio_sample_t) * nframes);
	    
	    //apply transformation
	    for (int i = 0; i < (int)nframes; i++) {
    		if (*(soft_clipper + i) >= 0) {
    			sgn = 1;
    		}
    		else {
    			sgn = -1;
    		}
    		
    		if((soft_clipper_gain * (*(soft_clipper + i))) < 0){
    		    temp = (-1) * (soft_clipper_gain * (*(soft_clipper + i)));
    		}
    		else{
    		    temp = (soft_clipper_gain * (*(soft_clipper + i)));
    		}
    		*(soft_clipper + i) = sgn * (1-(pow(e,(-1*temp))));
    	}
    	
    	out = jack_port_get_buffer (output_port, nframes);
		memcpy (out, soft_clipper, sizeof (jack_default_audio_sample_t) * nframes);
		free(soft_clipper);
    	
	}
	else if(strcmp(currentFilter, "Fuzz") == 0) {
	    
	    // variables set externally
        double e = 2.7182818;
        double temp = 0;
	    
	    // params
	    double fuzz_gain = 1.4; // the gain of the fuzz (range: > 1 to unsure)
	    
	    // allocate and init arrays
	    jack_default_audio_sample_t *samples;
		samples = (jack_default_audio_sample_t *)malloc(sizeof(jack_default_audio_sample_t ) * nframes);
		memcpy (samples, in, sizeof (jack_default_audio_sample_t) * nframes);
		
		jack_default_audio_sample_t *fuzz;
		fuzz = (jack_default_audio_sample_t *)malloc(sizeof(jack_default_audio_sample_t ) * nframes);
		for (int i = 0; i < (int)nframes; i++){
		    *(fuzz + i ) = 0;
		}
	    
	    jack_default_audio_sample_t *upsampled;
		upsampled = (jack_default_audio_sample_t *)malloc(sizeof(jack_default_audio_sample_t ) * (nframes * 8));
		for (int i = 0; i < ((int)nframes * 8); i++){
		    *(upsampled + i ) = 0;
		}
		
		jack_default_audio_sample_t *lowPass;
		lowPass = (jack_default_audio_sample_t *)malloc(sizeof(jack_default_audio_sample_t ) * nframes * 8);
		for (int i = 0; i < ((int)nframes * 8); i++){
		    *(lowPass + i) = 0;
		}
		
		jack_default_audio_sample_t *output;
		output = (jack_default_audio_sample_t *)malloc(sizeof(jack_default_audio_sample_t ) * nframes);
		for (int i = 0; i < (int)nframes; i++){
		    *(output + i) = 0;
		}
	    
	    // fuzzing transformation
	    for (int i = 0; i < (int)nframes; i++) {
	        if((*(fuzz + i)) < 0){
	            temp = (-1) * (*(samples + i));
	        }
	        else {
	            temp = (*(samples + i));
	        }
    		if(temp < 0.000001){
    		    *(fuzz + i) = 0.0;
    		}
    		else {
    		*(fuzz + i) = (*(samples + i) / temp) * (1 - pow(e, ((fuzz_gain * (pow((*(samples + i)), 2))) / temp)));
    		}
    	}
    	
    	// oversampling to reduce aliased frequencies
    	int oversampling = 0;
    	for(int i = 0; i < ((int)nframes * 8); i++){
    	    if(i % 8 == 0){
    	        *(upsampled + i) = *(fuzz + oversampling);
    	        oversampling++;
    	    }
    	    else{
    	        *(upsampled + i) = 0.0;
    	    }
    	} 
    	
    	// low pass transformation to take out high frequencies
    	// this value (d) should not be touched
    	double d = 0.6;
    	for (int i = 0; i < ((int)nframes * 8); i++) {
    		if (i == 0) {
    			*(lowPass + i) = *(upsampled + i);
    		}
    		else {
    			*(lowPass + i) = ((1.0 - d) * (*(upsampled + i))) + (d * (*(lowPass + i - 1)));
    		}
    	}
    	
    	// downsampling and sending to output array
    	int count = 0;
    	for(int i = 0; i < ((int)nframes * 8); i++){
    	    if(i % 8 ==0) {
    	        *(output + count) = *(lowPass + i);
    	        count++;
    	    }
    	}
    	
    	out = jack_port_get_buffer (output_port, nframes);
		memcpy (out, output, sizeof (jack_default_audio_sample_t) * nframes);
	    free(samples);
		free(fuzz);
		free(upsampled);
		free(lowPass);
		free(output);
	}
	else if(strcmp(currentFilter, "Overdrive") == 0) {
	    
	    jack_default_audio_sample_t *overdrive;
		overdrive = (jack_default_audio_sample_t *)malloc(sizeof(jack_default_audio_sample_t ) * nframes);
		memcpy (overdrive, in, sizeof (jack_default_audio_sample_t) * nframes);
	    
    	for (int i = 0; i < (int)nframes; i++) {
    		if (*(overdrive + i) > max_overdrive) {
    			*(overdrive + i) = max_overdrive;
    		}
    	}
    	
    	// no adjustable params here
    	
    	// overdrive transformation
    	for (int i = 0; i < (int)nframes; i++) {
    		if (*(overdrive + i) < ((1 / 3) * max_overdrive)) {
    			*(overdrive + i) = *(overdrive + i) * 2;
    		}
    		else if (*(overdrive + i) < ((2 / 3) * max_overdrive)) {
    			*(overdrive + i) = 1-(pow((2-(3*(*(overdrive + i)))),2)/3);
    		}
    		else {
    			*(overdrive + i) = *(overdrive + i);
    		}
    	}
    	
    	out = jack_port_get_buffer (output_port, nframes);
		memcpy (out, overdrive, sizeof (jack_default_audio_sample_t) * nframes);
		free(overdrive);
	}
	else if(strcmp(currentFilter, "Delay") == 0) {

		jack_default_audio_sample_t *delay;
		delay = (jack_default_audio_sample_t *)malloc(sizeof(jack_default_audio_sample_t ) * nframes);
		memcpy (delay, in, sizeof (jack_default_audio_sample_t) * nframes);

		int numSamples = (int)nframes;      // How many audio samples to process
      	int dpr = 0;                        // read pointer for delay buffer
    	int dpw = 0;                        // write pointer for delay buffer
    	if(initDelay == 0){
    	    for (int i = 0; i < delayBufLength; i++) {
            	*(delayData + i) = 0;
            	// if this doesnt work try *(delayData + i) = *(delay + i)
            }
            initDelay = 1;
    	}
    						
    	// User-adjustable effect parameters:
    	float dryMix_ = 0.5;      // Level of the dry (undelayed) signal (range: 0.1 to 0.9)
    	float wetMix_ = 0.5;      // Level of the wet (delayed) signal   (range: 0.1 to 0.9)
    	float feedback_ = 0.0;    // Feedback level (0 if no feedback)   (range: unsure)
    	for (int i = 0; i < numSamples; ++i)
    	{
    		const float in = *(delay + i);
    		float out = 0.0;
    		// The output is the input plus the contents of the 
    		// delay buffer (weighted by the mix levels).
    		out = (dryMix_ * in + wetMix_ * (*(delayData + dpr)));
    		// Store the current information in the delay buffer. 
    		// delayData[dpr] is the delay sample we just read, i.e. 
    		// what came out of the buffer. delayData[dpw] is what 
    		// we write to the buffer, i.e. what goes in
    		*(delayData + dpw) = in + (*(delayData + dpr) * feedback_);
    		if (++dpr >= delayBufLength) {
    			dpr = 0;
    		}
    		if (++dpw >= delayBufLength) { //increment by one then check
    			dpw = 0;
    		}
    		// Store output sample in buffer, replacing the input
    		*(delay + i) = out;  //*(samples + i) = out;
    	}
    		out = jack_port_get_buffer (output_port, nframes);
    		memcpy (out, delay, sizeof (jack_default_audio_sample_t) * nframes);
    		free(delay);
    	}
	else {
		printf("ERROR: Unkown Filter(%s)\n", currentFilter);
        out = jack_port_get_buffer (output_port, nframes);
		memcpy (out, in, sizeof (jack_default_audio_sample_t) * nframes);
	}

	return 0;      
}

/**
 * JACK calls this shutdown_callback if the server ever shuts down or
 * decides to disconnect the client.
 */
void
jack_shutdown (void *arg)
{
    // freeing circular buffers on shutdown
    free(chorusData);
    free(delayData);
    free(flangerData);
	exit (1);
}

int
main (int argc, char *argv[])
{
	
	const char **ports;
	const char *client_name = "simple";
	const char *server_name = NULL;
	jack_options_t options = JackNullOption;
	jack_status_t status;
	
	//initializing circular buffers outside of the process function so they keep their data 
    chorusData = (float*)malloc(sizeof(float) * chorusBufLength);
    delayData = (double*)malloc(sizeof(double) * delayBufLength);
    flangerData = (float*)malloc(sizeof(float) * flangerBufLength);
	
	/* open a client connection to the JACK server */
    currentFilter = "No Filter";
	client = jack_client_open (client_name, options, &status, server_name);
	if (client == NULL) {
		fprintf (stderr, "jack_client_open() failed, "
			 "status = 0x%2.0x\n", status);
		if (status & JackServerFailed) {
			fprintf (stderr, "Unable to connect to JACK server\n");
		}
		exit (1);
	}
	if (status & JackServerStarted) {
		fprintf (stderr, "JACK server started\n");
	}
	if (status & JackNameNotUnique) {
		client_name = jack_get_client_name(client);
		fprintf (stderr, "unique name `%s' assigned\n", client_name);
	}

	/* tell the JACK server to call `process()' whenever
	   there is work to be done.
	*/

	jack_set_process_callback (client, process, 0);

	/* tell the JACK server to call `jack_shutdown()' if
	   it ever shuts down, either entirely, or if it
	   just decides to stop calling us.
	*/

	jack_on_shutdown (client, jack_shutdown, 0);

	/* display the current sample rate. 
	 */

	printf ("engine sample rate: %" PRIu32 "\n",
		jack_get_sample_rate (client));

	/* create two ports */

	input_port = jack_port_register (client, "input",
					 JACK_DEFAULT_AUDIO_TYPE,
					 JackPortIsInput, 0);
	output_port = jack_port_register (client, "output",
					  JACK_DEFAULT_AUDIO_TYPE,
					  JackPortIsOutput, 0);

	if ((input_port == NULL) || (output_port == NULL)) {
		fprintf(stderr, "no more JACK ports available\n");
		exit (1);
	}

	/* Tell the JACK server that we are ready to roll.  Our
	 * process() callback will start running now. */

	if (jack_activate (client)) {
		fprintf (stderr, "cannot activate client");
		exit (1);
	}

	/* Connect the ports.  You can't do this before the client is
	 * activated, because we can't make connections to clients
	 * that aren't running.  Note the confusing (but necessary)
	 * orientation of the driver backend ports: playback ports are
	 * "input" to the backend, and capture ports are "output" from
	 * it.
	 */

	ports = jack_get_ports (client, NULL, NULL,
				JackPortIsPhysical|JackPortIsOutput);
	if (ports == NULL) {
		fprintf(stderr, "no physical capture ports\n");
		exit (1);
	}

	if (jack_connect (client, ports[0], jack_port_name (input_port))) {
		fprintf (stderr, "cannot connect input ports\n");
	}

	free (ports);
	
	ports = jack_get_ports (client, NULL, NULL,
				JackPortIsPhysical|JackPortIsInput);
	if (ports == NULL) {
		fprintf(stderr, "no physical playback ports\n");
		exit (1);
	}

	if (jack_connect (client, jack_port_name (output_port), ports[0])) {
		fprintf (stderr, "cannot connect output ports\n");
	}

	free (ports);

	/* keep running until stopped by the user */


	// Button listener
	/* Fire up threads */
    for (long i = 0; i < sizeof (threads) / sizeof (threadInfo_t); i++) {
        if (pthread_create(&threads[i].tid, NULL, service, (void *)i) < 0) {
            perror("pthread_create()");
            exit(-1);
        }
    }

    tcgetattr(0, &origtc);                         // get orig tty settings
    newtc = origtc;                                // copy them
    newtc.c_lflag &= ~ICANON;                      // put in '1 key mode'
    newtc.c_lflag &= ~ECHO;                        // turn off echo
	printf("Current Filter: %s\n", currentFilter);

    for(;;) {
        tcsetattr(0, TCSANOW, &newtc);             // echo off 1-key read mode
        char c = getchar();                        // get single key immed.
        tcsetattr(0, TCSANOW, &origtc);            // settings back to normal
        //printf("'%c'\n", c);        // show user what we got
        for (int i = 0; i < sizeof (threads) / sizeof (threadInfo_t); i++) {
            threadInfo_t *t = &threads[i];         // get shorthand ptr to thread
            if (strchr(t->keys, c) != NULL) {      // this thread listens for this key
                pthread_mutex_lock(&t->kqmutex);   // lock other threads out while we tamper 
                keyQueue_t *kq = calloc(sizeof (struct keyQueue), 1); // allocate pkt
                kq->key = c;                       // stash key there
                keyQueue_t *kptr = &t->kqhead;     // get pointer to queue head
                while(kptr->next != NULL)          // find first empty slot
                    kptr = kptr->next;
                kptr->next = kq;                   // enqueue key packet to thread
                pthread_mutex_unlock(&t->kqmutex); // unlock key queue
            }
        }
    }

	/* this is never reached but if the program
	   had some other way to exit besides being killed,
	   they would be important to call.
	*/

	jack_client_close (client);
	exit (0);
}
