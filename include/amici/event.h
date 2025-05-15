#ifndef AMICI_EVENT_H
#define AMICI_EVENT_H

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
     */
    Event(std::string id, bool initial_value)
        : id_(id)
        , initial_value_(initial_value) {}

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

/**
 * @brief Manages pending events.
 *
 * Stores events that trigger and are waiting to be handled.
 *
 * Until event priorities are implemented, this is a FIFO queue.
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

        // TODO: FIFO for now. Implement priority handling here.
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
