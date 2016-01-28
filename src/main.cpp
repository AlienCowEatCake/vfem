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

int paused(int ret)
{
#if defined _WIN32
    char * session = getenv("SESSIONNAME");
    if(!session || stricmp(session, "Console"))
        system("pause");
#endif
    return ret;
}

int main(int argc, char * argv [])
{
#if !defined _WIN32 && defined USE_NOSIGHUP
    // Устанавливаем обработчик SIGHUP
    struct sigaction sigact;
    memset(& sigact, 0, sizeof(struct sigaction));
    sigemptyset(& sigact.sa_mask);
    sigact.sa_handler = sighup_handler;
    sigaction(SIGHUP, & sigact, 0);
#endif

#if defined VFEM_USE_PML
    cout << " # VFEM_USE_PML" << endl;
#endif
#if defined USE_CXX11
    cout << " # USE_CXX11" << endl;
#endif
#if defined USE_NOSIGHUP
    cout << " # USE_NOSIGHUP" << endl;
#endif

    time_t seconds = time(NULL);
    char timebuf[24];
    strftime(timebuf, 24, "%Y-%m-%d_%H-%M-%S", localtime(&seconds));
    unsigned long time_exec = mtime();
    unsigned long time_solve = 0;

    bool nosolve = false;
    bool nopost = false;
    string config = "config.ini";
    string config_dir = "";
    for(int i = 1; i < argc; i++)
    {
        if(strcmp(argv[i], "-nosolve") == 0)        nosolve = true;
        else if(strcmp(argv[i], "-nopost") == 0)    nopost = true;
        else if(argv[i][0] != '-')
        {
            config = argv[i];
            size_t delim_pos = config.find_last_of("/");
#if defined _WIN32
            size_t delim_pos_w = config.find_last_of("\\");
            if(delim_pos == string::npos || (delim_pos_w != string::npos && delim_pos_w > delim_pos))
                delim_pos = delim_pos_w;
#endif
            if(delim_pos != string::npos)
                config_dir = config.substr(0, delim_pos + 1);
        }
        else cerr << "[Main] Unknown argument \"" << argv[i] << "\"" << endl;
    }

    VFEM v;
    if(!v.config.load(config)) return paused(1);
#if defined VFEM_USE_PML
    if(!v.config.load_pml(v.config.filename_pml)) return paused(1);
#endif
    if(!v.input_phys(config_dir + v.config.filename_phys)) return paused(1);
    if(!v.input_mesh(config_dir + v.config.filename_mesh)) return paused(1);
    v.make_struct();
    time_solve = mtime();
    if(!nosolve)
    {
        v.make_data();
        v.solve();
    }
    else
        if(!(v.config.filename_slae != "" && v.slae.restore_x(v.config.filename_slae)))
            return paused(1);
    time_solve = mtime() - time_solve;
    if(!nopost)
        postprocessing(v, timebuf);

    print_time(time_solve, "Solve time");
    time_exec = mtime() - time_exec;
    print_time(time_exec, "All time");
    return paused(0);
}
