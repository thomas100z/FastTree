from .node import Node


class AlignmentParser:
    def __init__(self, file_name: str) -> None:
        self.sequences = []

        try:
            with open(file_name) as file:
                lines = file.read().splitlines()
                for i in range(round(len(lines) / 2)):
                    self.sequences.append(Node(lines[(i * 2)].strip('>'), lines[(i * 2) + 1].strip()))
        except:
            raise Exception("Input file not in correct format!")

    def get_data(self) -> list:
        """
        Retrieves the alignment data in a Node object
        :return: list of Node
        """
        return self.sequences
