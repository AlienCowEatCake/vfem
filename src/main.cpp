#include "problems/problems.h"
#include <ctime>

#if defined _WIN32
#include <windows.h>
unsigned long mtime()
{
    return GetTickCount();
}
#else
unsigned long mtime()
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, & t);
    return (unsigned long)t.tv_sec * 1000 + t.tv_nsec / 1000000;
}
#endif

#if !defined _WIN32 && defined USE_NOSIGHUP
#include <sys/types.h>
#include <unistd.h>
#include <signal.h>

// Обработчик SIGHUP
void sighup_handler(int)
{
    // Обычно этот сигнал приходит, когда отключается терминал
    // Перехватим этот сигнал и перенаправим выходные потоки консоли
    // куда-нибудь в файлы, чтобы не потерять информацию или
    // проследить за процессом (не использовать при отладке!)
    char timebuf[24];
    time_t seconds = time(NULL);
    strftime(timebuf, 24, "%Y-%m-%d_%H-%M-%S", localtime(&seconds));
    stringstream ss;
    ss << getpid() << "_" << timebuf << ".txt";
    string fn_out = "stdout_" + ss.str();
    string fn_err = "stderr_" + ss.str();
    * stdout = * fopen(fn_out.c_str(), "w");
    * stderr = * fopen(fn_err.c_str(), "w");
    std::ios::sync_with_stdio();
}
#endif

void print_time(unsigned long msec, const string & descr)
{
    unsigned long seconds = msec / 1000;
    if(seconds > 3600)
    {
        unsigned long h = seconds / 3600;
        unsigned long m = (seconds - h * 3600) / 60;
        unsigned long s = seconds - h * 3600 - m * 60;
        cout << descr << ": \t" << h << " hr " << m << " min " << s << " sec." << endl;
    }
    else if(seconds > 60)
    {
        unsigned long m = seconds / 60;
        unsigned long s = seconds - m * 60;
        cout << descr << ": \t" << m << " min " << s << " sec." << endl;
    }
    else if(seconds > 0)
    {
        unsigned long ms = msec - seconds * 1000;
        cout << descr << ": \t" << seconds << " sec " << ms << " msec." << endl;
    }
    else
    {
        cout << descr << ": \t" << msec << " msec." << endl;
    }
}

int main()
{
#if !defined _WIN32 && defined USE_NOSIGHUP
    // Устанавливаем обработчик SIGHUP
    struct sigaction sigact;
    memset(& sigact, 0, sizeof(struct sigaction));
    sigemptyset(& sigact.sa_mask);
    sigact.sa_handler = sighup_handler;
    sigaction(SIGHUP, & sigact, 0);
#endif

    /**/
    cout << "Configuration:" << endl;
    cout << " # BASIS_ORDER: " << BASIS_ORDER << endl;
    cout << " # BASIS_TYPE:  " << BASIS_TYPE << endl;
#if defined VFEM_USE_PML
    cout << " # VFEM_USE_PML" << endl;
#endif
#if defined VFEM_USE_NONHOMOGENEOUS_FIRST
    cout << " # VFEM_USE_NONHOMOGENEOUS_FIRST" << endl;
#endif
#if defined VFEM_USE_ANALYTICAL
    cout << " # VFEM_USE_ANALYTICAL" << endl;
#endif
#if defined USE_CXX11
    cout << " # USE_CXX11" << endl;
#endif
#if defined USE_NOSIGHUP
    cout << " # USE_NOSIGHUP" << endl;
#endif
    cout << " # Mesh file: " << mesh_filename << endl;
    cout << " # Phys file: " << phys_filename << endl;
    /**/

    string new_tecplot_name = tecplot_filename.substr(0, tecplot_filename.find_last_of("."));
    time_t seconds = time(NULL);
    char timebuf[24];
    strftime(timebuf, 24, "%Y-%m-%d_%H-%M-%S", localtime(&seconds));
    new_tecplot_name += "_" + string(timebuf) + ".plt";
    unsigned long time_exec = mtime();
    unsigned long time_solve = 0;

    try
    {
        VFEM v;

        v.config.load("config.ini");
        system("pause");
        return 0;

        v.input_phys(phys_filename);
        v.input_mesh(mesh_filename);
        v.make();
        time_solve = mtime();
        v.solve();
        time_solve = mtime() - time_solve;
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
    time_exec = mtime() - time_exec;
    print_time(time_exec, "All time");
#if defined _WIN32
    system("pause");
#endif
    return 0;
}
