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
    /**
     * @brief Add a log entry
     * @param severity Severity level
     * @param identifier Short identifier for the logged event
     * @param message A more detailed message
     */
    void log(LogSeverity severity, std::string identifier, std::string message);

    /**
     * @brief Add a log entry with printf-like message formatting
     * @param severity Severity level
     * @param identifier Short identifier for the logged event
     * @param message A more detailed message
     */
    void log(LogSeverity severity, std::string identifier, const char* format, ...);

    /** The log items */
    std::vector<LogItem> items;
};


/**
 * @brief A log item.
 */
struct LogItem
{
    /**
     * @brief Construct a LogItem
     * @param severity
     * @param identifier
     * @param message
     */
    LogItem(
        LogSeverity severity,
        std::string identifier,
        std::string message
        ):
                  severity(severity)
                  ,identifier(identifier)
                  ,message(message)
                  {};

    /** Severity level */
    LogSeverity severity;

    /** Short identifier for the logged event */
    std::string identifier;

    /** A more detailed and readable message */
    std::string message;
};

} // namespace amici
#endif // AMICI_LOGGER_H
