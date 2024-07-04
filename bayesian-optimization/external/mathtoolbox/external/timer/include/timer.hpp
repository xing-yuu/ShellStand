/*
 MIT License
 
 Copyright (c) 2018 Yuki Koyama
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
 */

#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>

namespace timer
{
    class Timer
    {
    public:
        Timer(const std::string& message = "timer", bool show_destruction_message = true) :
        m_message(message),
        m_show_destruction_message(show_destruction_message),
        m_t_construct(std::chrono::system_clock::now())
        {
        }
        
        ~Timer()
        {
            if (m_show_destruction_message)
            {
                const std::string elapsed_time_message = get_elapsed_time_as_string();

                std::cout << "[ " << m_message << " : \t" << elapsed_time_message << " ]" << std::endl;
            }
        }

        std::string get_elapsed_time_as_string() const
        {
            const auto t_elapsed_in_microsec = get_elapsed_time_in_microseconds();
            const auto t_elapsed_in_millisec = t_elapsed_in_microsec / 1000.0;
            const auto t_elapsed_in_sec      = t_elapsed_in_millisec / 1000.0;

            std::ostringstream sstream;
            if (t_elapsed_in_sec > 10.0)
            {
                sstream << std::fixed << std::setprecision(3) << t_elapsed_in_sec << " s";
            }
            else if (t_elapsed_in_millisec > 10.0)
            {
                sstream << std::fixed << std::setprecision(3) << t_elapsed_in_millisec << " ms";
            }
            else
            {
                sstream << t_elapsed_in_microsec << " us";
            }
            return sstream.str();
        }

        long get_elapsed_time_in_microseconds() const
        {
            const auto t_now = std::chrono::system_clock::now();
            return std::chrono::duration_cast<std::chrono::microseconds>(t_now - m_t_construct).count();
        }

        long get_elapsed_time_in_milliseconds() const
        {
            const auto t_now = std::chrono::system_clock::now();
            return std::chrono::duration_cast<std::chrono::milliseconds>(t_now - m_t_construct).count();
        }

    private:
        const std::string                           m_message;
        const bool                                  m_show_destruction_message;
        const std::chrono::system_clock::time_point m_t_construct;
    };
}

#endif // TIMER_HPP
