#include "problems/problems.h"
#include <ctime>

#if !defined _WIN32 && defined USE_NOSIGHUP
#include <signal.h>

// Обработчик SIGHUP
void sighup_handler(int)
{
    // Обычно этот сигнал приходит, когда отключается терминал
    // Перехватим этот сигнал и перенаправим выходные потоки консоли
    // куда-нибудь в файлы, чтобы не потерять информацию или
    // проследить за процессом
    char timebuf[24];
    time_t seconds = time(NULL);
    strftime(timebuf, 24, "%Y-%m-%d_%H-%M-%S", localtime(&seconds));
    string fn_out = "stdout_" + string(timebuf) + ".txt";
    string fn_err = "stderr_" + string(timebuf) + ".txt";

    * stdout = * fopen(fn_out.c_str(), "w");
    * stderr = * fopen(fn_err.c_str(), "w");
    std::ios::sync_with_stdio();
}
#endif

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
#if defined VFEM_USE_PML_TENSOR
    cout << " # VFEM_USE_PML_TENSOR" << endl;
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
