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

// Распечатка затраченного времени
void print_time(unsigned long msec, const string & descr)
{
    cout << descr << ": \t" << msec << " msec \t(";
    unsigned long seconds = msec / 1000;
    if(seconds > 3600)
    {
        unsigned long h = seconds / 3600;
        unsigned long m = (seconds - h * 3600) / 60;
        unsigned long s = seconds - h * 3600 - m * 60;
        cout << h << " hr " << m << " min " << s << " sec";
    }
    else if(seconds > 60)
    {
        unsigned long m = seconds / 60;
        unsigned long s = seconds - m * 60;
        cout << m << " min " << s << " sec";
    }
    else if(seconds > 0)
    {
        unsigned long ms = msec - seconds * 1000;
        cout << seconds << " sec " << ms << " msec";
    }
    else
    {
        cout << msec << " msec";
    }
    cout << ")." << endl;
}

// Пауза, если это необходимо и поддерживаемо системой
int paused(int ret)
{
#if defined _WIN32
    char * session = getenv("SESSIONNAME");
    if(!session || stricmp(session, "Console"))
        system("pause");
#endif
    return ret;
}

// Решение одной задачи
bool vfem_solve(VFEM & v, const string & config, bool nosolve, bool nopost, const char * timebuf)
{
    unsigned long time_exec = mtime();
    unsigned long time_solve = 0;

    string config_dir = "";
    size_t delim_pos = config.find_last_of("/");
#if defined _WIN32
    size_t delim_pos_w = config.find_last_of("\\");
    if(delim_pos == string::npos || (delim_pos_w != string::npos && delim_pos_w > delim_pos))
        delim_pos = delim_pos_w;
#endif
    if(delim_pos != string::npos)
        config_dir = config.substr(0, delim_pos + 1);

    if(!v.config.load(config)) return false;
#if defined VFEM_USE_PML
    if(!v.config.load_pml(v.config.filename_pml)) return false;
#endif
    if(!v.input_phys(config_dir + v.config.filename_phys)) return false;
    if(!v.input_mesh(config_dir + v.config.filename_mesh)) return false;
    v.make_struct();
    if(!nosolve)
    {
        v.make_data();
        time_solve = mtime();
        v.solve();
        time_solve = mtime() - time_solve;
    }
    else
    {
        cout << "Restoring solution ..." << endl;
        if(v.config.filename_slae.empty())
        {
            cout << "Error in " << __FILE__ << ":" << __LINE__
                 << " empty VFEM::config.filename_slae" << endl;
            return false;
        }
        time_solve = mtime();
        if(!v.slae.restore_x(v.config.filename_slae))
            return false;
        time_solve = mtime() - time_solve;
    }
    if(!nopost)
        postprocessing(v, timebuf);

    print_time(time_solve, "Solve time");
    time_exec = mtime() - time_exec;
    print_time(time_exec, "All time");
    return true;
}

// Параметры командной строки для одной задачи
class cl_param
{
public:
    bool nosolve;
    bool nopost;
    string config;
    cl_param()
    {
        nosolve = false;
        nopost = false;
        config = "config.ini";
    }
};

// Main
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
#if defined USE_MKL
    cout << " # USE_MKL" << endl;
#endif

    time_t seconds = time(NULL);
    char timebuf[24];
    strftime(timebuf, 24, "%Y-%m-%d_%H-%M-%S", localtime(&seconds));

    vector<cl_param> cl_params;
    cl_param cl_param_curr;
    enum
    {
        DIFF_NO,
        DIFF_SIMPLE,
        DIFF_COMPLEX
    };
    int diff_type = DIFF_NO;

    for(int i = 1; i < argc; i++)
    {
        if(strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "-help") == 0 || strcmp(argv[i], "--help") == 0)
        {
            string name(argv[0]);
            size_t delim_pos = name.find_last_of("/");
#if defined _WIN32
            size_t delim_pos_w = name.find_last_of("\\");
            if(delim_pos == string::npos || (delim_pos_w != string::npos && delim_pos_w > delim_pos))
                delim_pos = delim_pos_w;
#endif
            if(delim_pos != string::npos)
                name = name.substr(delim_pos + 1);
            cout << name << " [-nosolve] [-nopost] [config.ini]" << endl;
            cout << name << " -diff [-nosolve] [-nopost] master.ini [-nosolve] [-nopost] slave.ini diff.ini" << endl;
            cout << name << " -diff_simple [-nosolve] [-nopost] master.ini [-nosolve] [-nopost] slave.ini diff.ini" << endl;
            cout << name << " -help" << endl;
            return paused(0);
        }
        else if(strcmp(argv[i], "-diff") == 0 || strcmp(argv[i], "--diff") == 0)
            diff_type = DIFF_COMPLEX;
        else if(strcmp(argv[i], "-diff_simple") == 0 || strcmp(argv[i], "--diff_simple") == 0)
            diff_type = DIFF_SIMPLE;
        else if(strcmp(argv[i], "-nosolve") == 0 || strcmp(argv[i], "--nosolve") == 0)
            cl_param_curr.nosolve = true;
        else if(strcmp(argv[i], "-nopost") == 0 || strcmp(argv[i], "--nopost") == 0)
            cl_param_curr.nopost = true;
        else if(argv[i][0] != '-')
        {
            cl_param_curr.config = argv[i];
            cl_params.push_back(cl_param_curr);
            cl_param_curr = cl_param();
        }
        else
        {
            cout << "[Main] Unknown argument \"" << argv[i] << "\"" << endl;
            return paused(1);
        }
    }

    if(diff_type != DIFF_NO && cl_params.size() != 3)
    {
        cout << "[Main] Incorrect command line arguments" << endl;
        return paused(1);
    }
    if(cl_params.size() == 0)
    {
        cl_params.push_back(cl_param_curr);
    }

    if(diff_type == DIFF_NO)
    {
        VFEM v;
        if(!vfem_solve(v, cl_params[0].config, cl_params[0].nosolve, cl_params[0].nopost, timebuf))
            return paused(1);
    }
    else
    {
        vector<diff_area> areas;
        if(!input_diff(cl_params[2].config, areas))
            return paused(1);
        VFEM master, slave;
        if(!vfem_solve(master, cl_params[0].config, cl_params[0].nosolve, cl_params[0].nopost, timebuf))
            return paused(1);
        if(!vfem_solve(slave, cl_params[1].config, cl_params[1].nosolve, cl_params[1].nopost, timebuf))
            return paused(1);
        if(diff_type == DIFF_SIMPLE)
            compare_simple(master, slave, areas);
        else
            compare_complex(master, slave, areas);
    }

    return paused(0);
}
