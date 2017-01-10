class Event:
    def __init__(self, dependent_addresses, impacted_addresses, priority=1):
        """
        A record and processor for a event type. Created by agent actions - contains method to update the grid with
        it's action
        :param dependent_addresses: Addresses required by this event
        :param impacted_addresses: Addresses this event updates
        :param priority: Can be used to order event occurrences
        """
        # Addresses which this event has a dependency on - changes to these addresses impacts whether the event happens
        self.dependent_addresses = dependent_addresses
        # Addresses which this event affects - their values will be amended in some way if this event happens
        self.impacted_addresses = impacted_addresses
        self.priority = priority

    def perform_event(self, work_grid):
        raise NotImplementedError
