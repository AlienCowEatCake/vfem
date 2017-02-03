#include "timers.h"
#include <ctime>
#include <iostream>

#if defined(_WIN32)
#include <windows.h>
#elif defined(__MACH__)
#include <mach/mach_time.h>
#endif

namespace fem_core { namespace utils { namespace timers {

/**
 * @brief Возвращает время в миллисекундах, само по себе ничего не значит, но по разнице можно засекать время
 * @return Время в миллисекундах
 */
unsigned long mtime()
{
#if defined(_WIN32)

    return GetTickCount();

#elif defined(__MACH__)

    mach_timebase_info_data_t timebase;
    mach_timebase_info(& timebase);
    uint64_t time = mach_absolute_time();
    return static_cast<unsigned long>
            ((time * static_cast<uint64_t>(timebase.numer)) /
            (static_cast<uint64_t>(timebase.denom) * static_cast<uint64_t>(1000000)));

#else

    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, & t);
    return static_cast<unsigned long>(t.tv_sec * 1000) +
           static_cast<unsigned long>(t.tv_nsec / 1000000);

#endif
}

/**
 * @brief Функция для распечатки затраченного времени, предназначена для работы совместно с mtime()
 * @param[in] msec Затраченное время в миллисекундах
 * @param[in] descr Текстовое описание операции
 */
void print_time(unsigned long msec, const std::string & descr)
{
    std::cout << descr << ": \t" << msec << " ms \t(";
    unsigned long seconds = msec / 1000;
    if(seconds > 604800)
    {
        unsigned long w = seconds / 604800;
        unsigned long d = (seconds - w * 604800) / 86400;
        unsigned long h = (seconds - w * 604800 - d * 86400) / 3600;
        unsigned long m = (seconds - w * 604800 - d * 86400 - h * 3600) / 60;
        unsigned long s = seconds - w * 604800 - d * 86400 - h * 3600 - m * 60;
        std::cout << w << " wk " << d << " day " << h << " hr " << m << " min " << s << " sec";
    }
    else if(seconds > 86400)
    {
        unsigned long d = seconds / 86400;
        unsigned long h = (seconds - d * 86400) / 3600;
        unsigned long m = (seconds - d * 86400 - h * 3600) / 60;
        unsigned long s = seconds - d * 86400 - h * 3600 - m * 60;
        std::cout << d << " day " << h << " hr " << m << " min " << s << " sec";
    }
    else if(seconds > 3600)
    {
        unsigned long h = seconds / 3600;
        unsigned long m = (seconds - h * 3600) / 60;
        unsigned long s = seconds - h * 3600 - m * 60;
        std::cout << h << " hr " << m << " min " << s << " sec";
    }
    else if(seconds > 60)
    {
        unsigned long m = seconds / 60;
        unsigned long s = seconds - m * 60;
        std::cout << m << " min " << s << " sec";
    }
    else if(seconds > 0)
    {
        unsigned long ms = msec - seconds * 1000;
        std::cout << seconds << " sec " << ms << " ms";
    }
    else
    {
        std::cout << msec << " ms";
    }
    std::cout << ")." << std::endl;
}

}}} // namespace fem_core::utils::timers
