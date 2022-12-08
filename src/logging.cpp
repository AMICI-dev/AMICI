#include "amici/logging.h"
#include "amici/misc.h"

#include <cstdarg>

namespace amici {

void Logger::log(LogSeverity severity, std::string identifier, std::string message)
{
    items.emplace_back(severity, identifier, message);
}

void Logger::log(LogSeverity severity, std::string identifier, const char *format, ...)
{
    va_list argptr;
    va_start(argptr, format);
    auto message = printfToString(format, argptr);
    va_end(argptr);

    log(severity, identifier, message);
}


} // namespace amici
