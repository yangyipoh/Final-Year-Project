//These codes were written for operation of MAXREFDES117 ppg sensor
//to collect raw heart rate signal

#include "Header.h"
#include <iterator>
// #include <dlib/statistics.h>

//using namespace arma;
using namespace std;

extern std::string directory;
extern std::ofstream rawheartrate; //Create a text file to store raw heart rate signal
extern std::ofstream rawheartrate_time;  //Create a text file to store time 
extern std::atomic<bool> isQuitRequested;

//To convert arma vec into stdvec (for tk::spline)
//typedef std::vector<double> stdvec;

//Raw heart rate signal collection function
void *hr_data_collection(void* threadarg) {
	//Create thread for this function
	struct hrthread_data *my_data;
    my_data = (struct hrthread_data *) threadarg;
	
	/*std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::chrono::milliseconds elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-my_data->startTime);*/
	
	//Record starting time
	rawheartrate_time << gpioTick() - my_data->startTime  << endl;
	
	while ( isQuitRequested == false)
	{	
		//Record data along with time
		rawheartrate << m.getIR() << endl;
		rawheartrate_time << gpioTick() - my_data->startTime << endl;
		/*end = std::chrono::steady_clock::now();
		elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end-my_data->startTime);*/
	}
	
	//Record ending time
	rawheartrate_time << gpioTick() - my_data->startTime << endl;

	rawheartrate.close();
	rawheartrate_time.close();
	pthread_exit(NULL);
}




