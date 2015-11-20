#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void readimg(char  *filename, float *img, int nx, int nz) {
	int size = nx * nz;
	FILE  *file = NULL;
	file = fopen(filename, "r");
	fread(img, sizeof(float), size, file);
	fclose(file);
}

void writeimg(char *filename, float *img, int nx, int nz) {
	int size = nx * nz;
	FILE  *file = NULL;
	file = fopen(filename, "w");
	fwrite(img, sizeof(float), size, file);
	fclose(file);
}

void fd2t10s_2d(int nx, int nz, const float* vel, float* prev_wave, float* curr_wave)
{
	float a[6];

	int d = 5;

	int ix, iz, k, curPos;

	a[5] = 1. / 3150.;
	a[4] = -5. / 1008.;
	a[3] = 5. / 126.;
	a[2] = -5. / 21.;
	a[1] = 5. / 3.;
	a[0] = a[1] + a[2] + a[3] + a[4] + a[5];

	for(ix = d; ix < nx - d; ix ++){
		for (iz = d; iz < nz - d; iz ++){
			curPos = ix * nz + iz;
			prev_wave[curPos] = 2 * curr_wave[curPos] - prev_wave[curPos];
			prev_wave[curPos] -= 4 * vel[curPos] * a[0] * curr_wave[curPos];
			for (k = 1; k <=5; k++){
				prev_wave[curPos] += vel[curPos] * a[k] * (curr_wave[curPos - k] + curr_wave[curPos + k] +
							curr_wave[curPos - k * nz] + curr_wave[curPos + k * nz]);
			}
		}
	}
}

void fd4t10s_2d(int nx, int nz, const float *vel, float *prev_wave, float *curr_wave) {
	float a[6];

	int d = 5;

	int ix, iz, curPos;

	a[0] = +1.53400796;
	a[1] = +1.78858721;
	a[2] = -0.31660756;
	a[3] = +0.07612173;
	a[4] = -0.01626042;
	a[5] = +0.00216736;

	float *u2 = malloc(sizeof(float) * nx * nz);

	for (ix = d; ix < nx - d; ix ++) {
		for (iz = d; iz < nz - d; iz ++) {
			curPos = ix * nz + iz;
			u2[curPos] = -4.0 * a[0] * curr_wave[curPos] +
				a[1] * (curr_wave[curPos - 1]  +  curr_wave[curPos + 1]  +
							curr_wave[curPos - nz]  +  curr_wave[curPos + nz])  +
				a[2] * (curr_wave[curPos - 2]  +  curr_wave[curPos + 2]  +
							curr_wave[curPos - 2 * nz]  +  curr_wave[curPos + 2 * nz])  +
				a[3] * (curr_wave[curPos - 3]  +  curr_wave[curPos + 3]  +
							curr_wave[curPos - 3 * nz]  +  curr_wave[curPos + 3 * nz])  +
				a[4] * (curr_wave[curPos - 4]  +  curr_wave[curPos + 4]  +
							curr_wave[curPos - 4 * nz]  +  curr_wave[curPos + 4 * nz])  +
				a[5] * (curr_wave[curPos - 5]  +  curr_wave[curPos + 5]  +
							curr_wave[curPos - 5 * nz]  +  curr_wave[curPos + 5 * nz]);
		}
	}

	for (ix = d + 1; ix < nx - d - 1; ix ++) { /// the range of ix is different from that in previous for loop
		for (iz = d + 1; iz < nz - d - 1; iz ++) { /// be careful of the range of iz
			curPos = ix * nz + iz;
			prev_wave[curPos] = 2.*curr_wave[curPos] - prev_wave[curPos] +
				(1.0f / vel[curPos]) * u2[curPos] + /// 2nd order
				1.0f / 12 * (1.0f / vel[curPos]) * (1.0f / vel[curPos]) *
				(u2[curPos - 1] + u2[curPos + 1] + u2[curPos - nz] + u2[curPos + nz] - 4 * u2[curPos]); /// 4th order
		}
	}
}

void forward_propagate_2d(int nx, int nz, float* vel, float* wav, int nt, float dt, int srcx, int srcz)
{

	int model_size = nx * nz;
	float* prev_wave = malloc(sizeof(float) * model_size);
	float* curr_wave = malloc(sizeof(float) * model_size);

	int ix;
	for(ix = 0; ix < model_size; ix ++)
	{
		prev_wave[ix] = 0.0f;
		curr_wave[ix] = 0.0f;
	}

	FILE* wave_fp = fopen("forward_wavefield.bin", "w");

	//propagation
	int time_step;
	for(time_step = 0; time_step < nt; time_step ++)
	{
		if (!(time_step%100)) printf("In Forward Propagate, time step: %d\n", time_step);
		// Inject Source
		int source_location = srcx * nz + srcz;
		curr_wave[source_location]  += wav[time_step];

		//One step propagation. And exchange the wave field buffer.
		fd2t10s_2d(nx, nz, vel, prev_wave, curr_wave);
		float*  tmp_wave = prev_wave;
		prev_wave = curr_wave;
		curr_wave = tmp_wave;

		if ((int)((time_step * dt)/(0.1)) != (int)((time_step-1)*dt/0.1))
		{
			fwrite(curr_wave, sizeof(float),  nx*nz, wave_fp);
		}

		int check_step = 5;

		if((time_step > 0) && (time_step != (nt - 1)) && !(time_step % check_step))
		{
			char check_file_name1[64];
			char check_file_name2[64];
			sprintf(check_file_name1, "/tmp/check_time_%d_1.su", time_step);
			sprintf(check_file_name2, "/tmp/check_time_%d_2.su", time_step);
			writeimg(check_file_name1, prev_wave, nx, nz);
			writeimg(check_file_name2, curr_wave, nx, nz);
		}
	}

	char check_file_name1[64] = "/tmp/check_time_last_1.su";
	char check_file_name2[64] = "/tmp/check_time_last_2.su";
	writeimg(check_file_name1, prev_wave, nx, nz);
	writeimg(check_file_name2, curr_wave, nx, nz);

	fclose(wave_fp);
	free(prev_wave);
	free(curr_wave);
}

void backward_propagate_onestep_2d(int nx, int nz, float* vel, float *wav, float* waveN, float* waveNm1, int nt, float dt, int srcx, int srcz, int time_step)
{
	int model_size = nx * nz;
	int check_step = 5;

	if (time_step  ==  nt - 1) {
		//Load last two time_step wave field
		char check_file_name1[64] = "/tmp/check_time_last_1.su";
		char check_file_name2[64] = "/tmp/check_time_last_2.su";
		readimg(check_file_name1, waveNm1, nx, nz);
		readimg(check_file_name2, waveN, nx, nz);
	}
	else if ((check_step > 0) && !(time_step % check_step) && (time_step != 0)) {
		char check_file_name1[64];
		char check_file_name2[64];
		sprintf(check_file_name1, "/tmp/check_time_%d_1.su", time_step);
		sprintf(check_file_name2, "/tmp/check_time_%d_2.su", time_step);
		readimg(check_file_name1, waveNm1, nx, nz);
		readimg(check_file_name2, waveN, nx, nz);
	}

	fd4t10s_2d(nx, nz, vel, waveN, waveNm1);
	float  *tmp = waveN;
	waveN = waveNm1;
	waveNm1 = tmp;

	if (time_step < nt) {
		int src_i = 0;
		for (src_i = 0; src_i < 1 ; src_i ++) {
			int source_location = srcx * nz + srcz;
			waveN[source_location] -= wav[time_step];
		}
	}
}

void receiver_propagate_onestep_2d(int nx, int nz, float* vel, float *wav, float* waveN, float* waveNm1, int nt, float dt, int srcx, int srcz, int time_step)
{
	if (time_step < nt) {
		int src_i = 0;
		for (src_i = 0; src_i < 1 ; src_i ++) {
			int source_location = srcx * nz + srcz;
			waveN[source_location] -= wav[time_step];
		}
	}

	fd2t10s_2d(nx, nz, vel, waveN, waveNm1);
	float  *tmp = waveN;
	waveN = waveNm1;
	waveNm1 = tmp;
}

int ricker_wavelet(float* wav, int nt, int st, float dt, float dx, float fpeak)
{

	int it;
	float x, xx;
	// A Ricker Wavelet
	for( it = 0; it < nt; it++)
	{
		x = 3.1415926536 * fpeak * ((it - st) * dt);
		xx = x * x;
		wav[it] = expf(-xx) * (1. - 2. * xx);
		wav[it] *= 1e6*dt*dt/dx/dx;
	}
	return 0;
}

int cross_correlation(float *src_wave, float *vsrc_wave, float *image, int model_size, float scale) {
	int i = 0;
	for (i = 0; i < model_size; i ++) {
		image[i] -= src_wave[i] * vsrc_wave[i] * scale;
	}
	return 0;
}

int main(int argc, char ** argv){
	// Init Parameters.
	int nx = 1000;
	int nz = 1000;
	float dx = 10.f;
	float dz = 10.f;
	int nt = 100;
	float freq = 10.f;
	float dt = 0.002f;
	int sx = 0;
	int sz = 0;
	float * vel;
	float * wav;

	// Init Velocity model and Source Wavelet.
	int model_size = nx * nz;
	// wave velocity
	vel = (float*) malloc(sizeof(float) * model_size);
	wav = (float*) malloc(sizeof(float) * nt);
	int ix, iz;
	for (ix = 0; ix < nx; ix ++){
		for (iz = 0; iz < nz; iz ++){
			vel[ix * nz + iz] = 2000. * 2000. * dt * dt / dx / dx;
		}
	}
	ricker_wavelet(wav, nt, 60, dt, dx, freq);
	FILE * wav_file = fopen("Wavelet.bin", "w");
	fwrite(wav, sizeof(float), nt, wav_file);

	// Do Forward Modeling
	forward_propagate_2d(nx, nz, vel, wav, nt, dt, sx, sz);

	float *N_src_wave = (float *)malloc(sizeof(float) * model_size);
	// N minus 1
	float *Nm1_src_wave = (float *)malloc(sizeof(float) * model_size);
	// virtual source
	float *N_vsrc_wave = (float *)malloc(sizeof(float) * model_size);
	float *Nm1_vsrc_wave = (float *)malloc(sizeof(float) * model_size);

	int idx;
	for (idx = 0; idx < model_size; idx++) {
		N_src_wave[idx] = 0.f;
		Nm1_src_wave[idx] = 0.f;
		N_vsrc_wave[idx] = 0.f;
		Nm1_vsrc_wave[idx] = 0.f;
	}

	int time_step;
	float *image;
	image = (float*) malloc(sizeof(float) * model_size);

	for (time_step = nt - 1; time_step >= 0; time_step--) {

		if (!(time_step % 100))
			printf("In Backward Propagate, time step: %d\n", time_step);
		backward_propagate_onestep_2d(nx, nz, vel, wav, N_src_wave, Nm1_src_wave, nt, dt, sx, sz, time_step);
		if (!(time_step % 100))
			printf("In Receiver Propagate, time step: %d\n", time_step);
		receiver_propagate_onestep_2d(nx, nz, vel, wav, N_vsrc_wave, Nm1_vsrc_wave, nt, dt, sx, sz, time_step);

		cross_correlation(N_src_wave, N_vsrc_wave, image, model_size, 1.0);
	}

	free(N_src_wave);
	free(Nm1_src_wave);
	free(N_vsrc_wave);
	free(Nm1_vsrc_wave);
	return 0;
}
