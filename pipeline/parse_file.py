from __future__ import annotations

class Profile:

    def __init__(self, sequences: ...) -> None:
        self.profile = Profile.form_profile_matrix_pseudocount(sequences)

    def get_profile(self) -> list[list[float]]:
        return self.profile

    @staticmethod
    def form_profile_matrix_pseudocount(matrix: list[str]) -> list[list[float]]:

        result = [[4 for j in range(len(matrix[0]))] for i in range(4)]

        base_amount = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        length = len(matrix)

        for i in range(len(matrix[0])):
            current_base_amount = dict(base_amount)

            for j in range(len(matrix)):
                current_base_amount[matrix[j][i]] += (1 / length)

            for k, score in enumerate(current_base_amount.values()):
                result[k][i] += score

        return result

    def join_profiles(self, other_profile: Profile) -> Profile:
        pass


class Alignment:

    def __init__(self, name: str, alignment: str) -> None:
        self.alignment = alignment
        self.name = name
        self.profile = Profile([alignment])


class AlignmentParser:

    def __init__(self, file_name: str) -> None:
        self.sequences = []
        lines = open(file_name).read().splitlines()

        for i in range(round(len(lines) / 2)):
            self.sequences.append(Alignment(lines[(i * 2)].strip(), lines[(i * 2) + 1].strip()))

    def get_data(self) -> list[Alignment]:
        return self.sequences
