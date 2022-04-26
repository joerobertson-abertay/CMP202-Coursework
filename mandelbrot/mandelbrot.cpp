// Mandelbrot set example
// Adam Sampson <a.sampson@abertay.ac.uk>

//#include <gl/GL.h>
//#include <gl/GLU.h>

#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <complex>
#include <fstream>
#include <iostream>
#include "amp.h"
#include <thread>
#include <mutex>
#include <condition_variable>
// Import things we need from the standard library
using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::complex;
using std::cout;
using std::endl;
using std::ofstream;
using std::thread;
using std::mutex;
using std::condition_variable;	
// Define the alias "the_clock" for the clock type we're going to use.
typedef std::chrono::steady_clock the_clock;

// The size of the image to generate.
const int WIDTH = 1920;
const int HEIGHT = 1200;
int sliceStart, sliceEnd = 0;
mutex coolMutex, tgaMutex;
condition_variable cv;
bool tgaRender = false;
// The number of times to iterate before we assume that a point isn't in the
// Mandelbrot set.
// (You may need to turn this up if you zoom further into the set.)
const int MAX_ITERATIONS = 500;

// The image data.
// Each pixel is represented as 0xRRGGBB.
uint32_t image[HEIGHT][WIDTH];


// Write the image to a TGA file with the given name.
// Format specification: http://www.gamers.org/dEngine/quake3/TGA.txt
void write_tga(const char *filename)
{
	std::unique_lock<mutex> tgaLock(tgaMutex);
	cv.wait(tgaLock, []() {return tgaRender; });
	ofstream outfile(filename, ofstream::binary);

	uint8_t header[18] = {
		0, // no image ID
		0, // no colour map
		2, // uncompressed 24-bit image
		0, 0, 0, 0, 0, // empty colour map specification
		0, 0, // X origin
		0, 0, // Y origin
		WIDTH & 0xFF, (WIDTH >> 8) & 0xFF, // width
		HEIGHT & 0xFF, (HEIGHT >> 8) & 0xFF, // height
		24, // bits per pixel
		0, // image descriptor
	};
	outfile.write((const char *)header, 18);

	for (int y = 0; y < HEIGHT; ++y)
	{
		for (int x = 0; x < WIDTH; ++x)
		{
			uint8_t pixel[3] = {
				image[y][x] & 0xFF, // blue channel
				(image[y][x] >> 8) & 0xFF, // green channel
				(image[y][x] >> 16) & 0xFF, // red channel
			};
			outfile.write((const char *)pixel, 3);
		}
	}

	outfile.close();
	if (!outfile)
	{
		// An error has occurred at some point since we opened the file.
		cout << "Error writing to " << filename << endl;
		exit(1);
	}
}


// Render the Mandelbrot set into the image array.
// The parameters specify the region on the complex plane to plot.
void compute_mandelbrot(double left, double right, double top, double bottom, int yStart, int yEnd)
{
	for (int y = yStart; y < yEnd; ++y)
	{
		for (int x = 0; x < WIDTH; ++x)
		{
			// Work out the point in the complex plane that
			// corresponds to this pixel in the output image.
			complex<double> c(left + (x * (right - left) / WIDTH),
				top + (y * (bottom - top) / HEIGHT));

			// Start off z at (0, 0).
			complex<double> z(0.0, 0.0);

			// Iterate z = z^2 + c until z moves more than 2 units
			// away from (0, 0), or we've iterated too many times.
			int iterations = 0;
			while (abs(z) < 2.0 && iterations < MAX_ITERATIONS)
			{
				z = (z * z) + c;

				++iterations;
			}

			if (iterations == MAX_ITERATIONS)
			{
				// z didn't escape from the circle.
				// This point is in the Mandelbrot set.
				image[y][x] = 0x000000; // black
			}
			else if (iterations >= 25)
			{
				image[y][x] = 0xBBBBBB; //grey
			}
			else
			{
				// z escaped within less than MAX_ITERATIONS
				// iterations. This point isn't in the set.
				image[y][x] = 0xFFFFFF; // white
			}
		}
	}
}

void tgaFunc()
{
	thread tgaThread(write_tga, "output.tga");
	tgaRender = true;
	tgaThread.join();
	cv.notify_one();
}

void mandelbrot_whatever()
{
	//locks for variables, just in case
	coolMutex.lock();
	sliceStart = sliceEnd;

	sliceEnd = sliceStart + 75;
	coolMutex.unlock();
	compute_mandelbrot(-2.0, 1.0, 1.0, -1.0, sliceStart, sliceEnd);
}

int main(int argc, char* argv[])
{
	cout << "Please wait..." << endl;

	// Start timing
	the_clock::time_point start = the_clock::now();

	std::vector<thread> threadVector;
	//const int threadCount = thread::hardware_concurrency();
	const int threadCount = 8; //this controls the number of threads
	// This shows the whole set.
	for (int i = 0; i < 16; i += threadCount) //16 strips
	{
		for (int j = 0; j < threadCount; j++) //makes threads
		{
			//important to use emplace_back over push_back here	
			threadVector.emplace_back(mandelbrot_whatever);
		}
		for (auto& k : threadVector)
		{
			k.join(); //joins threads
		};
		threadVector.clear(); //empties vector
	}

	// This zooms in on an interesting bit of detail.
	//compute_mandelbrot(-0.751085, -0.734975, 0.118378, 0.134488);

	// Stop timing
	the_clock::time_point end = the_clock::now();

	// Compute the difference between the two times in milliseconds
	auto time_taken = duration_cast<milliseconds>(end - start).count();
	cout << "Computing the Mandelbrot set took " << time_taken << " ms." << endl;
	
	tgaFunc();
	return 0;
}
