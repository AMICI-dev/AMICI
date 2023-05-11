#include "amici/logging.h"
#include "amici/misc.h"

#include <cstdarg>

namespace amici {

void Logger::log(
    LogSeverity severity, std::string const& identifier,
    std::string const& message
) {
    items.emplace_back(severity, identifier, message);
}

void Logger::log(
    LogSeverity severity, std::string const& identifier, char const* format, ...
) {
    va_list argptr;
    va_start(argptr, format);
    auto message = printfToString(format, argptr);
    va_end(argptr);

    log(severity, identifier, message);
}

} // namespace amici
