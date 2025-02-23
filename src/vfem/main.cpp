#include "common/common.h"
#include "problems/problems.h"

// Пауза, если это необходимо и поддерживаемо системой
int paused(int ret)
{
#if defined(_WIN32)
    char * session = getenv("SESSIONNAME");
    if(!session || compare_ci(session, "Console"))
        system("pause");
#endif
    return ret;
}

// Режимы решения
enum solve_mode
{
    MODE_NORMAL,
    MODE_NOSOLVE,
    MODE_CONTINUE
};

// Решение одной задачи
bool vfem_solve(VFEM & v, const string & config, solve_mode mode, bool nopost, const char * timebuf)
{
    cout << "Task configuration:" << endl;
    cout << " # mode = ";
    switch(mode)
    {
    case MODE_NORMAL:
        cout << "normal" << endl;
        break;
    case MODE_NOSOLVE:
        cout << "nosolve" << endl;
        break;
    case MODE_CONTINUE:
        cout << "continue" << endl;
        break;
    default:
        cout << "unknown" << endl;
        break;
    }
    cout << " # nopost = " << (nopost ? "true" : "false") << endl;

    unsigned long time_exec = mtime();
    unsigned long time_solve = 0;

    string config_dir = "";
    size_t delim_pos = config.find_last_of("/");
#if defined(_WIN32)
    size_t delim_pos_w = config.find_last_of("\\");
    if(delim_pos == string::npos || (delim_pos_w != string::npos && delim_pos_w > delim_pos))
        delim_pos = delim_pos_w;
#endif
    if(delim_pos != string::npos)
        config_dir = config.substr(0, delim_pos + 1);

    if(!v.config.load(config)) return false;
#if defined(VFEM_USE_PML)
    if(!v.config.load_pml(v.config.filename_pml)) return false;
#endif
    if(!v.input_phys(config_dir + v.config.filename_phys)) return false;
    if(!v.input_mesh(config_dir + v.config.filename_mesh)) return false;
    v.make_struct();
    if(mode == MODE_NORMAL)
    {
        v.make_data();
        time_solve = mtime();
        v.solve();
        time_solve = mtime() - time_solve;
    }
    else
    {
        if(mode == MODE_CONTINUE)
            v.make_data();
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
        if(mode == MODE_CONTINUE)
            v.solve();
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
    solve_mode mode;
    bool nopost;
    string config;
    cl_param()
    {
        mode = MODE_NORMAL;
        nopost = false;
        config = "config.ini";
    }
};

// Main
int main(int argc, char * argv [])
{
#if defined(USE_NOSIGHUP)
    set_nosighup();
#endif

    string build_conf;
#if defined(VFEM_USE_PML)
    build_conf.append(" # VFEM_USE_PML\n");
#endif
#if defined(USE_CXX11)
    build_conf.append(" # USE_CXX11\n");
#endif
#if defined(USE_NOSIGHUP)
    build_conf.append(" # USE_NOSIGHUP\n");
#endif
#if defined(USE_OMP)
    build_conf.append(" # USE_OMP\n");
#endif
#if defined(USE_MKL)
    build_conf.append(" # USE_MKL\n");
#endif
    if(!build_conf.empty())
        cout << "Build configuration:" << endl << build_conf << flush;

    string env_conf;
    char * env_curr;
    env_curr = getenv("NO_PROGRESS");
    if(env_curr)
    {
        env_conf.append(" # NO_PROGRESS = ");
        env_conf.append(env_curr);
        env_conf.append("\n");
    }
#if defined(USE_OMP)
    env_curr = getenv("OMP_NUM_THREADS");
    if(env_curr)
    {
        env_conf.append(" # OMP_NUM_THREADS = ");
        env_conf.append(env_curr);
        env_conf.append("\n");
    }
#endif
#if defined(USE_MKL)
    env_curr = getenv("MKL_NUM_THREADS");
    if(env_curr)
    {
        env_conf.append(" # MKL_NUM_THREADS = ");
        env_conf.append(env_curr);
        env_conf.append("\n");
    }
#endif
#if defined(_WIN32)
    env_curr = getenv("SESSIONNAME");
    if(env_curr)
    {
        env_conf.append(" # SESSIONNAME = ");
        env_conf.append(env_curr);
        env_conf.append("\n");
    }
#endif
    env_curr = getenv("GMRES_M");
    if(env_curr)
    {
        env_conf.append(" # GMRES_M = ");
        env_conf.append(env_curr);
        env_conf.append("\n");
    }
    if(!env_conf.empty())
        cout << "Environment configuration:" << endl << env_conf << flush;

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
#if defined(_WIN32)
            size_t delim_pos_w = name.find_last_of("\\");
            if(delim_pos == string::npos || (delim_pos_w != string::npos && delim_pos_w > delim_pos))
                delim_pos = delim_pos_w;
#endif
            if(delim_pos != string::npos)
                name = name.substr(delim_pos + 1);
            cout << name << " [-nosolve|-continue] [-nopost] [config.ini]" << endl;
            cout << name << " -diff [-nosolve|-continue] [-nopost] master.ini [-nosolve|-continue] [-nopost] slave.ini diff.ini" << endl;
            cout << name << " -diff_simple [-nosolve|-continue] [-nopost] master.ini [-nosolve|-continue] [-nopost] slave.ini diff.ini" << endl;
            cout << name << " -help" << endl;
            return paused(0);
        }
        else if(strcmp(argv[i], "-diff") == 0 || strcmp(argv[i], "--diff") == 0)
            diff_type = DIFF_COMPLEX;
        else if(strcmp(argv[i], "-diff_simple") == 0 || strcmp(argv[i], "--diff_simple") == 0)
            diff_type = DIFF_SIMPLE;
        else if(strcmp(argv[i], "-nosolve") == 0 || strcmp(argv[i], "--nosolve") == 0)
            cl_param_curr.mode = MODE_NOSOLVE;
        else if(strcmp(argv[i], "-continue") == 0 || strcmp(argv[i], "--continue") == 0)
            cl_param_curr.mode = MODE_CONTINUE;
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
        if(!vfem_solve(v, cl_params[0].config, cl_params[0].mode, cl_params[0].nopost, timebuf))
            return paused(1);
    }
    else
    {
        vector<diff_area> areas;
        if(!input_diff(cl_params[2].config, areas))
            return paused(1);
        VFEM master, slave;
        if(!vfem_solve(master, cl_params[0].config, cl_params[0].mode, cl_params[0].nopost, timebuf))
            return paused(1);
        if(!vfem_solve(slave, cl_params[1].config, cl_params[1].mode, cl_params[1].nopost, timebuf))
            return paused(1);
        if(diff_type == DIFF_SIMPLE)
            compare_simple(master, slave, areas);
        else
            compare_complex(master, slave, areas);
    }

    return paused(0);
}
