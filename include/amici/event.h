#ifndef AMICI_EVENT_H
#define AMICI_EVENT_H

#include "amici/defines.h"
#include <functional>
#include <list>
#include <stdexcept>
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
     * @param priority Priority of the event assignments
     */
    Event(std::string id, bool initial_value, realtype priority)
        : id_(id)
        , initial_value_(initial_value)
        , priority_(priority) {}

    /**
     * @brief Get the ID of the event
     * @return The ID of the event.
     */
    std::string const& get_id() const { return id_; }

    /**
     * @brief Get the initial value of the root function
     * @return The value of the root function at t_0.
     */
    bool get_initial_value() { return initial_value_; }

    /**
     * @brief Get the priority of the event assignments
     * @return The priority of the event assignments, or NAN if undefined.
     */
    realtype get_priority() const { return priority_; }

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

    /**
     * @brief The priority of the event assignments.
     *
     * When multiple events to be executed at the same time, their order is
     * determined by their priority. Higher priority event assignments are
     * executed first. If at least one priority is undefined (represented
     * by NaN), the order is undefined.
     */

    realtype priority_ = NAN;
};

/**
 * @brief Manages pending events.
 *
 * Stores events that trigger and are waiting to be handled.
 *
 * Manages the order of event execution based on their priority value.
 */
class EventQueue {
  public:
    /**
     * @brief EventQueue constructor
     */
    EventQueue() = default;

    /**
     * @brief Check if the queue is empty
     * @return True if the queue is empty, false otherwise
     */
    bool empty() const { return pending_events_.empty(); }

    /**
     * @brief Push an event to the queue
     * @param event The event to push
     */
    void push(Event const& event) {
        pending_events_.push_back(std::ref(event));
    }

    /**
     * @brief Get the next event to handle and remove it from the queue
     * @return The next event to handle
     * @throws std::runtime_error if there are no pending events
     */
    Event const& pop() {
        if (empty()) {
            throw std::runtime_error("No pending events");
        }

        // Sort events by priority from high to low.
        pending_events_.sort([](Event const& a, Event const& b) {
            // The priority is NaN if not defined. In this case, the execution
            // order is undefined, so this does not need any special treatment.
            return a.get_priority() > b.get_priority();
        });

        auto event = pending_events_.front();
        pending_events_.pop_front();
        return event.get();
    }

  private:
    /** The events waiting to be handled. */
    std::list<std::reference_wrapper<Event const>> pending_events_;
};
} // namespace amici

#endif // AMICI_EVENT_H
