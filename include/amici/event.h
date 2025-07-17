#ifndef AMICI_EVENT_H
#define AMICI_EVENT_H

#include "amici/defines.h"
#include "amici/model_state.h"

#include <algorithm>
#include <iterator>
#include <list>
#include <optional>
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
     * @param use_values_from_trigger_time Whether the event assignment is
     * evaluated using the state from the time point at which the event
     * triggered (true), or at the time point at which the event assignment
     * is evaluated (False).
     * @param initial_value Initial value of the root function
     * @param priority Priority of the event assignments
     */
    Event(
        std::string id, bool use_values_from_trigger_time, bool initial_value,
        realtype priority
    )
        : id_(id)
        , use_values_from_trigger_time_(use_values_from_trigger_time)
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
    bool get_initial_value() const { return initial_value_; }

    /**
     * @brief Get the priority of the event assignments
     * @return The priority of the event assignments, or NAN if undefined.
     */
    realtype get_priority() const { return priority_; }

    /**
     * @brief Check if the event assignment is evaluated using the state from
     * the time point at which the event is triggered or executed.
     * @return True if the event assignment is evaluated using the state from
     * the time point at which the event triggered, false otherwise.
     */
    bool uses_values_from_trigger_time() const {
        return use_values_from_trigger_time_;
    }

  private:
    /** The unique ID of this event. */
    std::string id_;

    /** Whether the event assignment is evaluated on the state from the time
     * point at which the event triggered (true), or at the time point at
     * which the event assignment is evaluated (false). */
    bool use_values_from_trigger_time_;

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
 * @brief PendingEvent struct
 *
 * Represents an event that has triggered and is waiting to be executed.
 */
struct PendingEvent {
    /** The event to be handled. */
    Event const& event;
    /** The index of `event` in the model */
    int idx;
    /**
     * The simulation state at the time `event` triggered.
     * Optional if `event.uses_values_from_trigger_time()` is false.
     */
    std::optional<SimulationState> state_old;
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
    void push(PendingEvent event) { pending_events_.push_back(event); }

    /**
     * @brief Get the next event to handle and remove it from the queue
     * @return The next event to handle
     * @throws std::runtime_error if there are no pending events
     */
    PendingEvent pop() {
        if (empty()) {
            throw std::runtime_error("No pending events");
        }

        // Sort events by priority from high to low.
        pending_events_.sort([](PendingEvent const& a, PendingEvent const& b) {
            // The priority is to NaN in AMICI if not defined.
            // In this case, the execution order is undefined,
            // so this does not need any special treatment.
            return a.event.get_priority() > b.event.get_priority();
        });

        // If there is any undefined (NaN) priority, the order of all
        // simulataneously executed events is undefined. We just proceed
        // with the first one in the list.
        if (std::any_of(
                pending_events_.begin(), pending_events_.end(),
                [](PendingEvent const& e) {
                    return std::isnan(e.event.get_priority());
                }
            )) {
            auto event = pending_events_.front();
            pending_events_.pop_front();
            return event;
        }

        // If there are multiple events with the same not-NaN priority,
        // then *their* ordering is random (not undefined).
        int num_equal_priority = 0;
        for (auto it = pending_events_.begin(); it != pending_events_.end();
             ++it) {
            if (it->event.get_priority()
                == pending_events_.front().event.get_priority()) {
                num_equal_priority++;
            } else {
                break;
            }
        }

        auto it = pending_events_.begin();
        if (num_equal_priority > 1) {
            // choose randomly among the events with the same priority
            std::advance(it, rand() % num_equal_priority);
        }

        auto event = *it;
        pending_events_.erase(it);
        return event;
    }

  private:
    /** The events waiting to be handled. */
    std::list<PendingEvent> pending_events_;
};
} // namespace amici

#endif // AMICI_EVENT_H
