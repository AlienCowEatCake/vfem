#include "nosighup.h"

#if !defined(_WIN32)

#include <sys/types.h>
#include <unistd.h>
#include <signal.h>
#include <ctime>
#include <sstream>
#include <string>
#include <cstring>
#include <cstdio>

namespace {

/**
 * @brief Обработчик SIGHUP
 */
void sighup_handler(int)
{
    // Обычно этот сигнал приходит, когда отключается терминал
    // Перехватим этот сигнал и перенаправим выходные потоки консоли
    // куда-нибудь в файлы, чтобы не потерять информацию или
    // проследить за процессом (не использовать при отладке!)
    char timebuf[24];
    time_t seconds = time(NULL);
    strftime(timebuf, 24, "%Y-%m-%d_%H-%M-%S", localtime(&seconds));
    std::stringstream ss;
    ss << getpid() << "_" << timebuf << ".txt";
    std::string fn_out = "stdout_" + ss.str();
    std::string fn_err = "stderr_" + ss.str();
    * stdout = * fopen(fn_out.c_str(), "w");
    * stderr = * fopen(fn_err.c_str(), "w");
    std::ios::sync_with_stdio();
}

} // namespace

#endif

namespace core { namespace utils { namespace nosighup {

/**
 * @brief Устанавливает обработчик SIGHUP с перенаправлением вывода в файлы
 * @note Hе использовать при отладке!
 */
void set_nosighup()
{
#if !defined(_WIN32)
    struct sigaction sigact;
    memset(& sigact, 0, sizeof(struct sigaction));
    sigemptyset(& sigact.sa_mask);
    sigact.sa_handler = sighup_handler;
    sigaction(SIGHUP, & sigact, 0);
#endif
}

}}} // namespace core::utils::nosighup
