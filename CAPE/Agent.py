class Agent:
    def __init__(self, address):
        """
        An Autonomous actor within the system. Abstract - should be subclassed.
        """
        self.address = address
        self.age = 0.0

    def output_code(self):
        """
        Code to be output for record grid. Should be overriden by subclass.
        :return:
        """
        raise NotImplementedError
