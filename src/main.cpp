#include "problems/problems.h"
#include <ctime>

void print_time(size_t seconds, const string & descr)
{
    if(seconds > 3600)
    {
        size_t h = (size_t)seconds / 3600;
        size_t m = (size_t)(seconds - h * 3600.0) / 60;
        size_t s = (size_t)(seconds - h * 3600.0 - m * 60.0);
        cout << descr << ": \t" << h << " hr " << m << " min " << s << " sec." << endl;
    }
    else if(seconds > 60)
    {
        size_t m = (size_t)(seconds) / 60;
        size_t s = (size_t)(seconds - m * 60.0);
        cout << descr << ": \t" << m << " min " << s << " sec." << endl;
    }
    else
    {
        size_t s = (size_t)seconds;
        cout << descr << ": \t" << s << " sec." << endl;
    }
}

int main()
{
    string new_tecplot_name = tecplot_filename.substr(0, tecplot_filename.find_last_of("."));
    time_t seconds = time(NULL);
    char timebuf[24];
    strftime(timebuf, 24, "%Y-%m-%d_%H-%M-%S", localtime(&seconds));
    new_tecplot_name += "_" + string(timebuf) + ".plt";
    size_t time_solve = 0;

    try
    {
        VFEM v;
        v.input_phys(phys_filename);
        v.input_mesh(mesh_filename);
        time_solve = (size_t)time(NULL);
        v.solve();
        time_solve = (size_t)time(NULL) - time_solve;
        v.output(new_tecplot_name);
        postprocessing(v, timebuf);
    }
    catch(int errn)
    {
        switch(errn)
        {
        case IO_FILE_ERROR:
            cerr << "Throw IO_FILE_ERROR" << endl;
            break;
        default:
            cerr << "Throw " << errn << " (unknown)" << endl;
            break;
        }
    }

    print_time(time_solve, "Solve time");
    size_t time_exec = (size_t)(time(NULL) - seconds);
    print_time(time_exec, "All time");
#if defined _WIN32
    system("pause");
#endif
    return 0;
}
