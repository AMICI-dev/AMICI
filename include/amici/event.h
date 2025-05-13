#ifndef AMICI_EVENT_H
#define AMICI_EVENT_H

#include <string>

namespace amici {

/**
 * @brief The Event class
 *
 * Represents an event. I.e., a potential discontinuity in the state or state
 * derivative that occurs when some trigger function transitions from `false` to
 * `true`.
 *
 * This will be extended to represent trigger time, trigger persistence,
 * event assignment priority, and other properties of the respective event
 * trigger and event assignments.
 */

class Event {
  public:
    /**
     * @brief Event constructor
     * @param id ID of the event
     * @param initial_value Initial value of the root function
     */
    Event(std::string id, bool initial_value)
        : id_(id)
        , initial_value_(initial_value) {}

    /**
     * @brief Get the initial value of the root function
     * @return The value of the root function at t_0.
     */
    bool get_initial_value() { return initial_value_; }

  private:
    /** The unique ID of this event. */
    std::string id_;

    /**
     * @brief Initial value of the trigger function.
     *
     * Events at t0 can only trigger if the initial value is set to `false`.
     * Must be specified during model compilation by setting the
     * `initialValue` attribute of an event trigger.
     */

    bool initial_value_ = true;
};

} // namespace amici

#endif // AMICI_EVENT_H
