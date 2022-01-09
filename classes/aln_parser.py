import numpy as np
from .node import Node


class AlignmentParser:
    def __init__(self, file_name: str) -> None:
        self.sequences = []

        with open(file_name) as file:
            lines = file.read().splitlines()
            for i in range(round(len(lines) / 2)):
                self.sequences.append(Node(lines[(i * 2)].strip('>'), lines[(i * 2) + 1].strip()))


    def get_data(self) -> list:
        return self.sequences
