#if defined(_MSC_VER) && !defined(_CRT_SECURE_NO_WARNINGS)
#define _CRT_SECURE_NO_WARNINGS
#endif
#include "progress.h"

#include <cstdlib>
#include <iostream>
#include <cstring>

namespace fem_core { namespace utils { namespace progress {

/**
 * @brief Распечатывает на stdout стильный прогрессбар
 * @param[in] message Заголовок прогрессбара (может быть пустым)
 * @param[in] index Индекс текущего обрабатываемого элемента
 * @param[in] num_all Общее количество обрабатываемых элементов
 */
void show_progress(const char * message, std::size_t index, std::size_t num_all)
{
    (void)(message);
    (void)(index);
    (void)(num_all);

    static const bool need_progress = std::getenv("NO_PROGRESS") == NULL;
    if(!need_progress) return;
/*
    const std::size_t tick = 1000;
    const char progress[4] = {'|', '/', '-', '\\'};
    index++;
    if(!(index % tick))
    {
        std::cout << "  ";
        if(std::strlen(message))
            std::cout << message << ": ";
        std::cout << index << progress[(index / tick) % 4] << num_all << '\r' << std::flush;
    }
    if(index >= num_all)
    {
        std::cout << "  ";
        if(std::strlen(message))
            std::cout << message << ": ";
        std::cout << num_all << progress[1] << num_all << std::endl;
    }
*/
/*
    static std::size_t perc_old = 100;
    static std::size_t diez_old = 72;

    if(index == 0 && std::strlen(message) > 0)
    {
        std::cout << "  * " << message << std::endl;
    }

    if(num_all > 1 && index < num_all - 1)
    {
        std::size_t norm = num_all - 1;
        std::size_t perc = (index * 100) / norm;
        std::size_t diez = (index * 72) / norm;

        if(perc_old != perc || diez_old != diez)
        {
            std::cout << '[';
            for(std::size_t i = 0; i < diez; i++)
                std::cout << '#';
            for(std::size_t i = diez; i < 72; i++)
                std::cout << ' ';
            std::cout << "] " << perc << "%\r" << std::flush;

            perc_old = perc;
            diez_old = diez;
        }
    }
    else
    {
        for(std::size_t i = 0; i < 79; i++)
            std::cout << ' ';
        std::cout << '\r' << std::flush;
    }
*/
/**/
    static std::size_t perc_old = 1000;
    static std::size_t psym_old = 0;
    const std::size_t syms_max = 79;
    const char progress[4] = {'|', '/', '-', '\\'};

    if(num_all > 1 && index < num_all - 1)
    {
        std::size_t norm = num_all - 1;
        //std::size_t perc = (index * 1000) / norm;
        //const std::size_t perc_len = 5;
        std::size_t perc = (index * 100) / norm;
        const std::size_t perc_len = 3;

        if(perc_old != perc)
        {
            std::size_t message_len = (message ? std::strlen(message) : 0);
            if(message_len)
                std::cout << message << ": |";
            else
                std::cout << '|';

            std::size_t syms_avail = syms_max - message_len - (message_len ? 3 : 1) - perc_len - 2;
            std::size_t equals = (index * syms_avail) / norm;
            for(std::size_t i = 0; i < equals; i++)
                std::cout << '=';
            for(std::size_t i = equals; i < syms_avail; i++)
                std::cout << ' ';
            //std::cout << progress[psym_old] << ' ' << perc / 10 << '.' << perc % 10 << "%\r" << std::flush;
            std::cout << progress[psym_old] << ' ' << perc << "%\r" << std::flush;

            perc_old = perc;
            if(psym_old == 3)
                psym_old = 0;
            else
                psym_old++;
        }
    }
    else
    {
        for(std::size_t i = 0; i < syms_max; i++)
            std::cout << ' ';
        std::cout << '\r' << std::flush;
    }
/**/
}

}}} // namespace fem_core::utils::progress
