#ifndef AMICI_LOGGER_H
#define AMICI_LOGGER_H

#include <string>
#include <vector>

namespace amici {

struct LogItem;

/**
 * @brief Severity levels for logging.
 */
enum class LogSeverity {
    error,
    warning,
    debug,
};

/**
 * @brief A logger, holding a list of error messages.
 */
class Logger
{
  public:
    Logger() = default;
    void log(LogSeverity severity, std::string identifier, std::string message);
    void log(LogSeverity severity, std::string identifier, const char* format, ...);

    std::vector<LogItem> items;
};


/**
 * @brief A log item.
 */
struct LogItem
{
    LogItem(
        LogSeverity severity,
        std::string identifier,
        std::string message
        ):
                  severity(severity)
                  ,identifier(identifier)
                  ,message(message)
                  {};
    LogSeverity severity;
    std::string identifier;
    std::string message;
};

} // namespace amici
#endif // AMICI_LOGGER_H
